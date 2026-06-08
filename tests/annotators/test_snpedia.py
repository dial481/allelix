# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Tests for the SNPedia annotator."""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.databases.snpedia_parser import _PARSER_VERSION

if TYPE_CHECKING:
    from pathlib import Path

from allelix.annotators.snpedia import SNPediaAnnotator
from allelix.models import Variant


@pytest.fixture()
def snpedia_data_dir(tmp_path: Path) -> Path:
    """Create a minimal SNPedia structured database for testing."""
    db_path = tmp_path / "snpedia.sqlite"
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        conn.executescript("""
            CREATE TABLE snpedia_genotypes (
                rsid TEXT NOT NULL,
                allele1 TEXT NOT NULL,
                allele2 TEXT NOT NULL,
                magnitude REAL,
                repute TEXT,
                summary TEXT,
                gene TEXT,
                scraped_at TEXT
            );
            CREATE INDEX idx_snpedia_rsid_alleles
                ON snpedia_genotypes(rsid, allele1, allele2);
            CREATE TABLE database_versions (
                name TEXT PRIMARY KEY,
                source_url TEXT NOT NULL,
                version TEXT,
                downloaded_at TEXT NOT NULL,
                record_count INTEGER NOT NULL,
                remote_signal TEXT
            );
        """)
        genotypes = [
            (
                "rs1801133",
                "C",
                "C",
                0.0,
                "Good",
                "Common genotype: normal homocysteine levels",
                "MTHFR",
                "2026-05-20T00:00:00",
            ),
            (
                "rs1801133",
                "C",
                "T",
                2.2,
                "Bad",
                "1 copy of C677T allele of MTHFR",
                "MTHFR",
                "2026-05-20T00:00:00",
            ),
            (
                "rs1801133",
                "T",
                "T",
                2.8,
                "Bad",
                "homozygous C677T of MTHFR",
                "MTHFR",
                "2026-05-20T00:00:00",
            ),
            (
                "rs4680",
                "A",
                "G",
                None,
                None,
                "Intermediate dopamine levels",
                "COMT",
                "2026-05-20T00:00:00",
            ),
            ("rs9999999", "A", "A", 0.0, "Good", None, None, "2026-05-20T00:00:00"),
            (
                "rs52820871",
                "G",
                "G",
                0.0,
                "Good",
                "common genotype",
                "TNFRSF13B",
                "2026-05-20T00:00:00",
            ),
            (
                "rs52820871",
                "G",
                "T",
                3.0,
                "Bad",
                "TACI variant",
                "TNFRSF13B",
                "2026-05-20T00:00:00",
            ),
            ("i3000001", "A", "G", 4.0, "Bad", "CF carrier", "CFTR", "2026-05-20T00:00:00"),
        ]
        conn.executemany(
            "INSERT INTO snpedia_genotypes "
            "(rsid, allele1, allele2, magnitude, repute, summary, gene, scraped_at) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            genotypes,
        )
        conn.execute(
            "INSERT INTO database_versions "
            "(name, source_url, version, downloaded_at, record_count, remote_signal) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (
                "snpedia",
                "https://bots.snpedia.com/api.php",
                "scraped 2026-05-20 (8 genotypes)",
                "2026-05-20T00:00:00",
                8,
                f"|pv:{_PARSER_VERSION}",
            ),
        )
        conn.commit()
    return tmp_path


class TestAnnotatorLifecycle:
    """Annotator setup, readiness, version, and close."""

    def test_is_ready(self, snpedia_data_dir: Path) -> None:
        ann = SNPediaAnnotator(snpedia_data_dir)
        assert ann.is_ready()
        ann.close()

    def test_not_ready_missing_db(self, tmp_path: Path) -> None:
        ann = SNPediaAnnotator(tmp_path)
        assert not ann.is_ready()

    def test_not_ready_no_structured_table(self, tmp_path: Path) -> None:
        db_path = tmp_path / "snpedia.sqlite"
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "CREATE TABLE pages (title TEXT PRIMARY KEY, category TEXT, "
                "content TEXT, scraped_at TEXT)"
            )
            conn.execute("INSERT INTO pages VALUES ('Rs1', 'snp', 'content', '2026-01-01')")
            conn.commit()
        ann = SNPediaAnnotator(tmp_path)
        assert not ann.is_ready()
        ann.close()

    def test_version(self, snpedia_data_dir: Path) -> None:
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = ann.version()
        assert v is not None
        assert "scraped" in v
        assert "8 genotypes" in v
        ann.close()

    def test_record_count(self, snpedia_data_dir: Path) -> None:
        ann = SNPediaAnnotator(snpedia_data_dir)
        count = ann.record_count()
        assert count == 8
        ann.close()

    def test_setup_is_noop(self, snpedia_data_dir: Path) -> None:
        ann = SNPediaAnnotator(snpedia_data_dir)
        ann.setup()
        ann.close()

    def test_requires_download_false(self) -> None:
        assert SNPediaAnnotator.requires_download is False

    def test_remote_signals_none(self, snpedia_data_dir: Path) -> None:
        ann = SNPediaAnnotator(snpedia_data_dir)
        assert ann.fetch_remote_signal() is None
        assert ann.cached_remote_signal() is None
        ann.close()


class TestAnnotateGenotype:
    """Genotype matching and annotation output."""

    def test_heterozygous_match(self, snpedia_data_dir: Path) -> None:
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(
            rsid="rs1801133",
            chromosome="1",
            position=11796321,
            allele1="C",
            allele2="T",
        )
        results = ann.annotate(v)
        assert len(results) == 1
        assert results[0].source == "snpedia"
        assert results[0].magnitude == 2.2
        assert results[0].gene == "MTHFR"
        assert results[0].attribution == "SNPedia"
        assert results[0].significance == "snpedia_bad"
        assert results[0].category == "clinical"
        ann.close()

    def test_homozygous_match(self, snpedia_data_dir: Path) -> None:
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(
            rsid="rs1801133",
            chromosome="1",
            position=11796321,
            allele1="T",
            allele2="T",
        )
        results = ann.annotate(v)
        assert len(results) == 1
        assert results[0].magnitude == 2.8
        ann.close()

    def test_reference_homozygous(self, snpedia_data_dir: Path) -> None:
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(
            rsid="rs1801133",
            chromosome="1",
            position=11796321,
            allele1="C",
            allele2="C",
        )
        results = ann.annotate(v)
        assert len(results) == 1
        assert results[0].magnitude == 0.0
        assert results[0].significance == "snpedia_good"
        assert results[0].category == "trait"
        ann.close()

    def test_allele_order_independent(self, snpedia_data_dir: Path) -> None:
        """User alleles T/C should match stored C/T row."""
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(
            rsid="rs1801133",
            chromosome="1",
            position=11796321,
            allele1="T",
            allele2="C",
        )
        results = ann.annotate(v)
        assert len(results) == 1
        assert results[0].magnitude == 2.2
        ann.close()

    def test_no_call_returns_empty(self, snpedia_data_dir: Path) -> None:
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(
            rsid="rs1801133",
            chromosome="1",
            position=11796321,
            allele1="-",
            allele2="-",
        )
        assert ann.annotate(v) == []
        ann.close()

    def test_unknown_rsid_returns_empty(self, snpedia_data_dir: Path) -> None:
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(rsid="rs999999999", chromosome="1", position=1, allele1="A", allele2="A")
        assert ann.annotate(v) == []
        ann.close()

    def test_genotype_mismatch_returns_empty(self, snpedia_data_dir: Path) -> None:
        """User has A/A but SNPedia only has C/C, C/T, T/T for rs1801133."""
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(
            rsid="rs1801133",
            chromosome="1",
            position=11796321,
            allele1="A",
            allele2="A",
        )
        assert ann.annotate(v) == []
        ann.close()

    def test_i_probe_annotation(self, snpedia_data_dir: Path) -> None:
        """I-prefixed 23andMe probes match against SNPedia I-probe data."""
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(rsid="i3000001", chromosome="7", position=117199646, allele1="A", allele2="G")
        results = ann.annotate(v)
        assert len(results) == 1
        assert results[0].description == "SNPedia: CF carrier"
        assert results[0].gene == "CFTR"
        assert any("snpedia.com" in r and "I3000001" in r for r in results[0].references)
        ann.close()

    def test_i_probe_wrong_alleles_returns_empty(self, snpedia_data_dir: Path) -> None:
        """I-probe with non-matching alleles produces no annotation."""
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(rsid="i3000001", chromosome="7", position=117199646, allele1="T", allele2="T")
        assert ann.annotate(v) == []
        ann.close()

    def test_i_probe_not_in_database_returns_empty(self, snpedia_data_dir: Path) -> None:
        """I-probe rsid not present in the database at all."""
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(rsid="i9999999", chromosome="1", position=1, allele1="A", allele2="G")
        assert ann.annotate(v) == []
        ann.close()

    def test_i_probe_uppercase_input_matches(self, snpedia_data_dir: Path) -> None:
        """Uppercase I3000001 from user input normalizes to lowercase for lookup."""
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(rsid="I3000001", chromosome="7", position=117199646, allele1="A", allele2="G")
        results = ann.annotate(v)
        assert len(results) == 1
        assert results[0].gene == "CFTR"
        ann.close()

    def test_non_rs_non_i_rsid_returns_empty(self, snpedia_data_dir: Path) -> None:
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(rsid="x12345", chromosome="1", position=1, allele1="A", allele2="A")
        assert ann.annotate(v) == []
        ann.close()

    def test_missing_magnitude_defaults_zero(self, snpedia_data_dir: Path) -> None:
        """rs4680(A;G) has no magnitude — should default to 0.0."""
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(rsid="rs4680", chromosome="22", position=19963748, allele1="A", allele2="G")
        results = ann.annotate(v)
        assert len(results) == 1
        assert results[0].magnitude == 0.0
        assert results[0].gene == "COMT"
        ann.close()

    def test_empty_summary_skipped(self, snpedia_data_dir: Path) -> None:
        """rs9999999(A;A) has no summary — should be skipped."""
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(rsid="rs9999999", chromosome="1", position=1, allele1="A", allele2="A")
        assert ann.annotate(v) == []
        ann.close()

    def test_description_attributed(self, snpedia_data_dir: Path) -> None:
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(
            rsid="rs1801133",
            chromosome="1",
            position=11796321,
            allele1="C",
            allele2="T",
        )
        results = ann.annotate(v)
        assert results[0].description.startswith("SNPedia:")
        ann.close()

    def test_references_include_snpedia_url(self, snpedia_data_dir: Path) -> None:
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(
            rsid="rs1801133",
            chromosome="1",
            position=11796321,
            allele1="C",
            allele2="T",
        )
        results = ann.annotate(v)
        assert any("snpedia.com" in r for r in results[0].references)
        ann.close()

    def test_gene_populated(self, snpedia_data_dir: Path) -> None:
        """rs52820871 should have gene TNFRSF13B from structured data."""
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(rsid="rs52820871", chromosome="13", position=1, allele1="G", allele2="T")
        results = ann.annotate(v)
        assert len(results) == 1
        assert results[0].magnitude == 3.0
        assert results[0].description == "SNPedia: TACI variant"
        assert results[0].gene == "TNFRSF13B"
        ann.close()

    def test_genotype_match_format(self, snpedia_data_dir: Path) -> None:
        """genotype_match should be concatenated alleles (e.g., 'CT')."""
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(
            rsid="rs1801133",
            chromosome="1",
            position=11796321,
            allele1="C",
            allele2="T",
        )
        results = ann.annotate(v)
        assert results[0].genotype_match == "CT"
        ann.close()


class TestADR0023HomRefSuppression:
    """ADR-0023: suppress SNPedia disease claims when user is hom-ref per ClinVar."""

    def test_suppresses_hom_ref(self, snpedia_data_dir: Path) -> None:
        """User hom-ref per ClinVar must not get SNPedia disease annotation."""
        ann = SNPediaAnnotator(
            snpedia_data_dir,
            clinvar_ref_provider=lambda rsid, build: "C" if rsid == "rs1801133" else None,
        )
        v = Variant(rsid="rs1801133", chromosome="1", position=11796321, allele1="C", allele2="C")
        assert ann.annotate(v) == []
        ann.close()

    def test_emits_when_user_carries_alt(self, snpedia_data_dir: Path) -> None:
        """User carrying ALT allele still gets the annotation."""
        ann = SNPediaAnnotator(
            snpedia_data_dir,
            clinvar_ref_provider=lambda rsid, build: "C" if rsid == "rs1801133" else None,
        )
        v = Variant(rsid="rs1801133", chromosome="1", position=11796321, allele1="C", allele2="T")
        results = ann.annotate(v)
        assert len(results) == 1
        assert results[0].magnitude == 2.2
        ann.close()

    def test_emits_when_no_clinvar_ref(self, snpedia_data_dir: Path) -> None:
        """rsid not in ClinVar — no suppression, SNPedia annotation emits."""
        ann = SNPediaAnnotator(
            snpedia_data_dir,
            clinvar_ref_provider=lambda rsid, build: None,
        )
        v = Variant(rsid="rs1801133", chromosome="1", position=11796321, allele1="C", allele2="T")
        results = ann.annotate(v)
        assert len(results) == 1
        ann.close()

    def test_emits_when_provider_is_none(self, snpedia_data_dir: Path) -> None:
        """No provider wired — behaves like pre-ADR-0023 (no suppression)."""
        ann = SNPediaAnnotator(snpedia_data_dir, clinvar_ref_provider=None)
        v = Variant(rsid="rs1801133", chromosome="1", position=11796321, allele1="C", allele2="C")
        results = ann.annotate(v)
        assert len(results) == 1
        ann.close()

    def test_hom_ref_suppression_het_passes(self, snpedia_data_dir: Path) -> None:
        """Het user is not hom-ref — annotation emits even with REF provider."""
        ann = SNPediaAnnotator(
            snpedia_data_dir,
            clinvar_ref_provider=lambda rsid, build: "T" if rsid == "rs1801133" else None,
        )
        v = Variant(rsid="rs1801133", chromosome="1", position=11796321, allele1="C", allele2="T")
        results = ann.annotate(v)
        assert len(results) == 1
        ann.close()


class TestSummarySuppressionFilter:
    """Suppress SNPedia annotations whose summary flags orientation uncertainty."""

    def test_mis_oriented_summary_suppressed(self, snpedia_data_dir: Path) -> None:
        """Pages whose summary flags orientation uncertainty must not emit."""
        db_path = snpedia_data_dir / "snpedia.sqlite"
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "INSERT INTO snpedia_genotypes VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    "rs1064651",
                    "C",
                    "C",
                    8.0,
                    "Bad",
                    "Gaucher disease, but more likely a mis-oriented interpretation",
                    "GBA",
                    "2026-05-20T00:00:00",
                ),
            )
            conn.commit()
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(rsid="rs1064651", chromosome="1", position=155204239, allele1="C", allele2="C")
        assert ann.annotate(v) == []
        ann.close()

    def test_normal_summary_not_suppressed(self, snpedia_data_dir: Path) -> None:
        """Normal summaries without suppression keywords emit normally."""
        ann = SNPediaAnnotator(snpedia_data_dir)
        v = Variant(rsid="rs1801133", chromosome="1", position=11796321, allele1="C", allele2="T")
        results = ann.annotate(v)
        assert len(results) == 1
        ann.close()


class TestAnnotatorGracefulAbsence:
    """Annotator degrades gracefully when SNPedia data is not available."""

    def test_not_ready_no_crash(self, tmp_path: Path) -> None:
        ann = SNPediaAnnotator(tmp_path)
        assert not ann.is_ready()
        assert ann.version() is None
        assert ann.record_count() is None
        ann.close()
