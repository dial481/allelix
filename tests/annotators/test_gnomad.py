# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Tests for the gnomAD population frequency annotator."""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.gnomad import _BULK_BATCH_SIZE, GnomadAnnotator
from allelix.databases.gnomad_loader import GNOMAD_DB_FILENAME
from allelix.databases.schema import GNOMAD_SCHEMA
from allelix.models import Annotation, Variant

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture
def gnomad_db(tmp_path: Path) -> Path:
    """Create a minimal gnomAD SQLite cache with test data."""
    db_path = tmp_path / GNOMAD_DB_FILENAME
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        for stmt in GNOMAD_SCHEMA.split(";"):
            stmt = stmt.strip()
            if stmt:
                conn.execute(stmt)
        conn.executemany(
            "INSERT INTO gnomad_frequencies"
            " (chrom, pos, ref, alt, rsid, af) VALUES (?, ?, ?, ?, ?, ?)",
            [
                ("1", 11796321, "G", "A", "rs1801133", 0.35),
                ("22", 19963748, "G", "A", "rs4680", 0.50),
                ("19", 44908684, "T", "C", "rs429358", 0.15),
                ("19", 44908684, "T", "G", "rs429358", 0.08),
                ("17", 43093517, "A", "G", "rs80357906", 0.0005),
                ("17", 43093610, "C", "T", "rs199476100", 0.000005),
            ],
        )
        conn.execute(
            "INSERT INTO database_versions"
            " (name, source_url, version, downloaded_at, record_count)"
            " VALUES (?, ?, ?, ?, ?)",
            ("gnomad", "test://mock", "4.1", "2026-01-01T00:00:00Z", 5),
        )
        conn.commit()
    return tmp_path


@pytest.fixture
def gnomad_full_db(tmp_path: Path) -> Path:
    """Create a gnomAD cache with >2M records (simulated via record_count)."""
    db_path = tmp_path / GNOMAD_DB_FILENAME
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        for stmt in GNOMAD_SCHEMA.split(";"):
            stmt = stmt.strip()
            if stmt:
                conn.execute(stmt)
        conn.execute(
            "INSERT INTO database_versions"
            " (name, source_url, version, downloaded_at, record_count)"
            " VALUES (?, ?, ?, ?, ?)",
            ("gnomad", "test://mock-full", "4.1", "2026-01-01T00:00:00Z", 16_000_000),
        )
        conn.commit()
    return tmp_path


class TestSetupAndStatus:
    """Annotator lifecycle: ready, version, record_count, close."""

    def test_unconfigured_is_not_ready(self, tmp_path: Path) -> None:
        annotator = GnomadAnnotator(tmp_path)
        assert not annotator.is_ready()

    def test_configured_is_ready(self, gnomad_db: Path) -> None:
        annotator = GnomadAnnotator(gnomad_db)
        assert annotator.is_ready()

    def test_version_returns_string(self, gnomad_db: Path) -> None:
        annotator = GnomadAnnotator(gnomad_db)
        assert annotator.version() == "4.1"

    def test_version_returns_none_when_missing(self, tmp_path: Path) -> None:
        annotator = GnomadAnnotator(tmp_path)
        assert annotator.version() is None

    def test_record_count(self, gnomad_db: Path) -> None:
        annotator = GnomadAnnotator(gnomad_db)
        assert annotator.record_count() == 5

    def test_record_count_returns_none_when_missing(self, tmp_path: Path) -> None:
        annotator = GnomadAnnotator(tmp_path)
        assert annotator.record_count() is None


class TestAnnotateReturnsEmpty:
    """gnomAD does not participate in the per-variant annotation loop."""

    def test_annotate_returns_empty_list(self, gnomad_db: Path) -> None:
        annotator = GnomadAnnotator(gnomad_db)
        v = Variant(
            rsid="rs1801133",
            chromosome="1",
            position=11796321,
            allele1="G",
            allele2="A",
            build="GRCh38",
        )
        assert annotator.annotate(v) == []

    def test_annotate_never_produces_annotations(self, gnomad_db: Path) -> None:
        """Even for known rsIDs, annotate() always returns []."""
        annotator = GnomadAnnotator(gnomad_db)
        for rsid in ["rs1801133", "rs4680", "rs429358"]:
            v = Variant(rsid=rsid, chromosome="1", position=1, allele1="A", allele2="A")
            assert annotator.annotate(v) == []


class TestLookup:
    """Single-rsID frequency lookup."""

    def test_known_rsid(self, gnomad_db: Path) -> None:
        annotator = GnomadAnnotator(gnomad_db)
        try:
            assert annotator.lookup("rs1801133") == pytest.approx(0.35)
        finally:
            annotator.close()

    def test_unknown_rsid(self, gnomad_db: Path) -> None:
        annotator = GnomadAnnotator(gnomad_db)
        try:
            assert annotator.lookup("rs999999999") is None
        finally:
            annotator.close()

    def test_rare_variant(self, gnomad_db: Path) -> None:
        annotator = GnomadAnnotator(gnomad_db)
        try:
            assert annotator.lookup("rs199476100") == pytest.approx(0.000005)
        finally:
            annotator.close()


class TestMultiAllelic:
    """MAX(af) aggregation for multi-allelic sites sharing an rsID."""

    def test_lookup_returns_max_af(self, gnomad_db: Path) -> None:
        """rs429358 has two alleles (af=0.15 and af=0.08); lookup returns 0.15."""
        annotator = GnomadAnnotator(gnomad_db)
        try:
            assert annotator.lookup("rs429358") == pytest.approx(0.15)
        finally:
            annotator.close()

    def test_bulk_lookup_returns_max_af(self, gnomad_db: Path) -> None:
        """bulk_lookup returns the MAX(af) across alleles for each rsID."""
        annotator = GnomadAnnotator(gnomad_db)
        try:
            result = annotator.bulk_lookup({"rs429358", "rs4680"})
            assert result["rs429358"] == pytest.approx(0.15)
            assert result["rs4680"] == pytest.approx(0.50)
        finally:
            annotator.close()


class TestBulkLookup:
    """Batched frequency lookup for annotation enrichment."""

    def test_all_found(self, gnomad_db: Path) -> None:
        annotator = GnomadAnnotator(gnomad_db)
        try:
            result = annotator.bulk_lookup({"rs1801133", "rs4680"})
            assert result == {"rs1801133": pytest.approx(0.35), "rs4680": pytest.approx(0.50)}
        finally:
            annotator.close()

    def test_partial_match(self, gnomad_db: Path) -> None:
        annotator = GnomadAnnotator(gnomad_db)
        try:
            result = annotator.bulk_lookup({"rs1801133", "rs999999999"})
            assert "rs1801133" in result
            assert "rs999999999" not in result
        finally:
            annotator.close()

    def test_empty_input(self, gnomad_db: Path) -> None:
        annotator = GnomadAnnotator(gnomad_db)
        try:
            assert annotator.bulk_lookup(set()) == {}
        finally:
            annotator.close()

    def test_no_matches(self, gnomad_db: Path) -> None:
        annotator = GnomadAnnotator(gnomad_db)
        try:
            result = annotator.bulk_lookup({"rs111111111", "rs222222222"})
            assert result == {}
        finally:
            annotator.close()

    def test_batching_over_limit(self, gnomad_db: Path) -> None:
        """Bulk lookup with more rsIDs than _BULK_BATCH_SIZE still works."""
        annotator = GnomadAnnotator(gnomad_db)
        try:
            rsids = {f"rs{i}" for i in range(1, _BULK_BATCH_SIZE + 100)}
            rsids.add("rs1801133")
            result = annotator.bulk_lookup(rsids)
            assert "rs1801133" in result
        finally:
            annotator.close()

    def test_bulk_batch_size_constant(self) -> None:
        """Pin: batch size stays within SQLite's variable limit."""
        assert _BULK_BATCH_SIZE == 900


class TestCloseable:
    """Context manager and close semantics."""

    def test_close_is_idempotent(self, gnomad_db: Path) -> None:
        annotator = GnomadAnnotator(gnomad_db)
        annotator.lookup("rs1801133")
        annotator.close()
        annotator.close()

    def test_context_manager(self, gnomad_db: Path) -> None:
        with GnomadAnnotator(gnomad_db) as annotator:
            assert annotator.is_ready()

    def test_context_manager_closes_connection(self, gnomad_db: Path) -> None:
        with GnomadAnnotator(gnomad_db) as annotator:
            annotator.lookup("rs1801133")
        assert annotator._conn is None


class TestRemoteSignal:
    """Remote signal fetch and cached signal retrieval."""

    def test_cached_signal_none_when_missing(self, tmp_path: Path) -> None:
        annotator = GnomadAnnotator(tmp_path)
        assert annotator.cached_remote_signal() is None

    def test_cached_signal_none_when_no_signal_stored(self, gnomad_db: Path) -> None:
        annotator = GnomadAnnotator(gnomad_db)
        assert annotator.cached_remote_signal() is None

    def test_cached_signal_returns_stored_value(self, gnomad_db: Path) -> None:
        db_path = gnomad_db / GNOMAD_DB_FILENAME
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET remote_signal = ? WHERE name = 'gnomad'",
                ("etag:abc123",),
            )
            conn.commit()
        annotator = GnomadAnnotator(gnomad_db)
        assert annotator.cached_remote_signal() == "etag:abc123"


class TestRegistryMetadata:
    """Class attributes match the annotator contract."""

    def test_name(self) -> None:
        assert GnomadAnnotator.name == "gnomad"

    def test_display_name(self) -> None:
        assert GnomadAnnotator.display_name == "gnomAD"

    def test_attribution(self) -> None:
        assert GnomadAnnotator.attribution == "gnomAD"

    def test_requires_download(self) -> None:
        assert GnomadAnnotator.requires_download is True


class TestPipelineEnrichment:
    """gnomAD frequency stamps on annotations via the pipeline."""

    def test_allele_frequency_stamped(self, gnomad_db: Path) -> None:
        """bulk_lookup results appear on annotation.allele_frequency."""
        annotator = GnomadAnnotator(gnomad_db)
        freq_map = annotator.bulk_lookup({"rs1801133", "rs4680", "rs999999999"})
        annotations = [
            Annotation(
                source="clinvar",
                rsid="rs1801133",
                significance="clinvar_pathogenic",
                category="clinical",
                magnitude=9.0,
                description="test",
                attribution="ClinVar",
                genotype_match="A/A",
            ),
            Annotation(
                source="clinvar",
                rsid="rs999999999",
                significance="clinvar_vus",
                category="clinical",
                magnitude=5.0,
                description="unknown rsid",
                attribution="ClinVar",
                genotype_match="T/T",
            ),
        ]
        for a in annotations:
            a.allele_frequency = freq_map.get(a.rsid)
        assert annotations[0].allele_frequency == pytest.approx(0.35)
        assert annotations[1].allele_frequency is None
        annotator.close()


class TestInstallPrebuiltCache:
    """install_prebuilt_cache decompresses and stamps signal."""

    def test_decompress_and_stamp(self, tmp_path: Path) -> None:
        import gzip

        from allelix.databases.gnomad_loader import install_prebuilt_cache

        src_db = tmp_path / "source.sqlite"
        with contextlib.closing(sqlite3.connect(src_db)) as conn:
            for stmt in GNOMAD_SCHEMA.split(";"):
                stmt = stmt.strip()
                if stmt:
                    conn.execute(stmt)
            conn.execute(
                "INSERT INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count)"
                " VALUES (?, ?, ?, ?, ?)",
                ("gnomad", "test://prebuilt", "4.1", "2026-01-01T00:00:00Z", 100),
            )
            conn.commit()

        gz_path = tmp_path / "test.sqlite.gz"
        with src_db.open("rb") as f_in, gzip.open(gz_path, "wb") as f_out:
            f_out.write(f_in.read())

        dest_db = tmp_path / "dest" / GNOMAD_DB_FILENAME
        dest_db.parent.mkdir()
        install_prebuilt_cache(gz_path, dest_db, remote_signal="etag:test123")

        assert dest_db.exists()
        with contextlib.closing(sqlite3.connect(dest_db)) as conn:
            row = conn.execute(
                "SELECT remote_signal FROM database_versions WHERE name = 'gnomad'"
            ).fetchone()
        assert row[0] == "etag:test123"

    def test_decompress_without_signal(self, tmp_path: Path) -> None:
        import gzip

        from allelix.databases.gnomad_loader import install_prebuilt_cache

        src_db = tmp_path / "source.sqlite"
        with contextlib.closing(sqlite3.connect(src_db)) as conn:
            for stmt in GNOMAD_SCHEMA.split(";"):
                stmt = stmt.strip()
                if stmt:
                    conn.execute(stmt)
            conn.execute(
                "INSERT INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count)"
                " VALUES (?, ?, ?, ?, ?)",
                ("gnomad", "test://prebuilt", "4.1", "2026-01-01T00:00:00Z", 100),
            )
            conn.commit()

        gz_path = tmp_path / "test.sqlite.gz"
        with src_db.open("rb") as f_in, gzip.open(gz_path, "wb") as f_out:
            f_out.write(f_in.read())

        dest_db = tmp_path / GNOMAD_DB_FILENAME
        install_prebuilt_cache(gz_path, dest_db)

        assert dest_db.exists()
        with contextlib.closing(sqlite3.connect(dest_db)) as conn:
            row = conn.execute(
                "SELECT remote_signal FROM database_versions WHERE name = 'gnomad'"
            ).fetchone()
        assert row[0] is None
