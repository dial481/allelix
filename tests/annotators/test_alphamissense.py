# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Tests for the AlphaMissense enrichment annotator."""

from __future__ import annotations

import sqlite3
from typing import TYPE_CHECKING

from allelix.annotators.alphamissense import AlphaMissenseAnnotator
from allelix.databases.schema import ALPHAMISSENSE_SCHEMA
from allelix.models import Variant

if TYPE_CHECKING:
    from pathlib import Path


def _build_db(tmp_path: Path) -> Path:
    """Build a minimal AlphaMissense SQLite cache for testing."""
    db_path = tmp_path / "alphamissense.sqlite"
    conn = sqlite3.connect(db_path)
    conn.executescript(ALPHAMISSENSE_SCHEMA)
    conn.executemany(
        "INSERT INTO alphamissense_scores "
        "(chrom, pos, ref, alt, rsid, uniprot_id, transcript_id, "
        "protein_variant, am_pathogenicity, am_class) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [
            (
                "1",
                100,
                "A",
                "G",
                "rs1001",
                "P12345",
                "ENST001",
                "A100G",
                0.95,
                "likely_pathogenic",
            ),
            ("1", 200, "C", "T", "rs1002", "P12345", "ENST001", "C200T", 0.20, "likely_benign"),
            ("2", 300, "G", "A", "rs2001", "P67890", "ENST002", "G300A", 0.45, "ambiguous"),
            ("3", 400, "T", "C", None, "P99999", "ENST003", "T400C", 0.80, "likely_pathogenic"),
        ],
    )
    conn.execute(
        "INSERT INTO database_versions (name, source_url, version, downloaded_at, record_count) "
        "VALUES ('alphamissense', 'https://test', '2023.1', '2026-06-08', 4)",
    )
    conn.commit()
    conn.close()
    return db_path


class TestAlphaMissenseAnnotator:
    def test_is_ready(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        assert am.is_ready()

    def test_not_ready_without_db(self, tmp_path: Path) -> None:
        am = AlphaMissenseAnnotator(tmp_path)
        assert not am.is_ready()

    def test_version(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        assert am.version() == "2023.1"

    def test_record_count(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        assert am.record_count() == 4

    def test_annotate_returns_empty(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        v = Variant(rsid="rs1001", chromosome="1", position=100, allele1="A", allele2="G")
        assert am.annotate(v) == []

    def test_lookup_found(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        result = am.lookup("rs1001")
        assert result is not None
        score, cls = result
        assert score == 0.95
        assert cls == "likely_pathogenic"
        am.close()

    def test_lookup_not_found(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        assert am.lookup("rs9999") is None
        am.close()

    def test_bulk_lookup(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        result = am.bulk_lookup({"rs1001", "rs1002", "rs9999"})
        assert len(result) == 2
        assert result["rs1001"] == (0.95, "likely_pathogenic")
        assert result["rs1002"] == (0.20, "likely_benign")
        assert "rs9999" not in result
        am.close()

    def test_bulk_lookup_empty(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        assert am.bulk_lookup(set()) == {}
        am.close()

    def test_close(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        am.lookup("rs1001")
        am.close()
        assert am._conn is None

    def test_am_class_thresholds(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        result = am.bulk_lookup({"rs1001", "rs1002", "rs2001"})
        assert result["rs1001"][1] == "likely_pathogenic"
        assert result["rs1002"][1] == "likely_benign"
        assert result["rs2001"][1] == "ambiguous"
        am.close()


class TestMultiAllelicMax:
    """MAX(am_pathogenicity) aggregation for multi-allelic sites sharing an rsID."""

    def test_bulk_lookup_returns_max_score(self, tmp_path: Path) -> None:
        db_path = tmp_path / "alphamissense.sqlite"
        conn = sqlite3.connect(db_path)
        conn.executescript(ALPHAMISSENSE_SCHEMA)
        conn.executemany(
            "INSERT INTO alphamissense_scores "
            "(chrom, pos, ref, alt, rsid, uniprot_id, transcript_id, "
            "protein_variant, am_pathogenicity, am_class) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                ("1", 100, "A", "G", "rs5000", "P1", "E1", "A1G", 0.20, "likely_benign"),
                ("1", 100, "A", "T", "rs5000", "P1", "E1", "A1T", 0.92, "likely_pathogenic"),
                ("1", 100, "A", "C", "rs5000", "P1", "E1", "A1C", 0.45, "ambiguous"),
            ],
        )
        conn.execute(
            "INSERT INTO database_versions "
            "(name, source_url, version, downloaded_at, record_count) "
            "VALUES ('alphamissense', 'test', '1.0', '2026-01-01', 3)",
        )
        conn.commit()
        conn.close()

        am = AlphaMissenseAnnotator(tmp_path)
        result = am.bulk_lookup({"rs5000"})
        assert result["rs5000"][0] == 0.92
        am.close()

    def test_lookup_returns_max_score(self, tmp_path: Path) -> None:
        db_path = tmp_path / "alphamissense.sqlite"
        conn = sqlite3.connect(db_path)
        conn.executescript(ALPHAMISSENSE_SCHEMA)
        conn.executemany(
            "INSERT INTO alphamissense_scores "
            "(chrom, pos, ref, alt, rsid, uniprot_id, transcript_id, "
            "protein_variant, am_pathogenicity, am_class) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                ("2", 200, "C", "T", "rs6000", "P2", "E2", "C2T", 0.10, "likely_benign"),
                ("2", 200, "C", "A", "rs6000", "P2", "E2", "C2A", 0.88, "likely_pathogenic"),
            ],
        )
        conn.execute(
            "INSERT INTO database_versions "
            "(name, source_url, version, downloaded_at, record_count) "
            "VALUES ('alphamissense', 'test', '1.0', '2026-01-01', 2)",
        )
        conn.commit()
        conn.close()

        am = AlphaMissenseAnnotator(tmp_path)
        result = am.lookup("rs6000")
        assert result is not None
        assert result[0] == 0.88
        am.close()
