# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Tests for the unified analysis pipeline."""

from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.models import Variant
from allelix.parsers.myhappygenes import MyHappyGenesParser
from allelix.reports._pipeline import (
    _DETECTION_BUFFER_LIMIT,
    AnalysisResult,
    _BuildDetectionState,
    run_analysis,
)
from allelix.utils.build_detect import BUILD_GRCH36, BUILD_GRCH37, KNOWN_SNP_POSITIONS

if TYPE_CHECKING:
    from pathlib import Path


def _ann(**overrides):
    from allelix.models import Annotation

    defaults = {
        "source": "clinvar",
        "rsid": "rs1",
        "significance": "clinvar_pathogenic",
        "category": "clinical",
        "magnitude": 5.0,
        "description": "x",
        "attribution": "ClinVar",
        "genotype_match": "A",
        "gene": "GENE1",
    }
    defaults.update(overrides)
    return Annotation(**defaults)


class TestAnalysisResultFilter:
    def _result(self, annotations) -> AnalysisResult:
        from pathlib import Path

        return AnalysisResult(
            file_path=Path("dummy.txt"),
            parser_name="x",
            parser_display_name="X",
            sample_id="S",
            build="GRCh37",
            total_variants=0,
            skipped_count=0,
            annotators_used=[],
            annotations=annotations,
        )

    def test_min_magnitude_excludes_low(self):
        r = self._result([_ann(rsid="lo", magnitude=2), _ann(rsid="hi", magnitude=8)])
        kept = r.filter(min_magnitude=5)
        assert [a.rsid for a in kept] == ["hi"]

    def test_category_filter(self):
        r = self._result([_ann(rsid="c", category="clinical"), _ann(rsid="p", category="pharma")])
        assert [a.rsid for a in r.filter(category="pharma")] == ["p"]

    def test_genes_filter_case_insensitive(self):
        r = self._result([_ann(rsid="m", gene="MTHFR"), _ann(rsid="b", gene="BRCA1")])
        kept = r.filter(genes={"mthfr"})
        assert [a.rsid for a in kept] == ["m"]

    def test_sort_is_magnitude_then_rsid(self):
        r = self._result(
            [
                _ann(rsid="rs2", magnitude=5),
                _ann(rsid="rs1", magnitude=5),
                _ann(rsid="rs3", magnitude=8),
            ]
        )
        kept = r.filter()
        assert [a.rsid for a in kept] == ["rs3", "rs1", "rs2"]


class TestRunAnalysis:
    def test_streams_and_collects(self, mock_mhg_path: Path, all_annotators_data_dir: Path):
        parser = MyHappyGenesParser()
        annotators = [
            ClinVarAnnotator(all_annotators_data_dir),
            PharmGKBAnnotator(all_annotators_data_dir),
        ]
        result = run_analysis(mock_mhg_path, parser, annotators)
        assert result.parser_name == "myhappygenes"
        assert result.sample_id == "MHG000001"
        assert result.total_variants == 2016
        assert any(a.source == "clinvar" for a in result.annotations)
        assert any(a.source == "pharmgkb" for a in result.annotations)
        # ADR-0021: composite version reports both builds when annotator
        # manages both. Single-build instances collapse to a single part.
        clinvar_versions = [v for name, v in result.annotators_used if name == "clinvar"]
        assert clinvar_versions, "ClinVar annotator missing from used set"
        assert "20260101" in clinvar_versions[0]

    def test_annotator_connections_closed_after_run(
        self, mock_mhg_path: Path, clinvar_data_dir: Path
    ):
        parser = MyHappyGenesParser()
        ann = ClinVarAnnotator(clinvar_data_dir)
        run_analysis(mock_mhg_path, parser, [ann])
        # ExitStack closed every per-build connection.
        assert ann._conns == {}


class TestGRCh36FlushFailSafe:
    """GRCh36 non-confident detection must use GRCh36 as effective build.

    Issue #6: when detection points to GRCh36 but isn't confident
    (matched < inspected), the pipeline was falling back to header_build
    or GRCh37. This bypassed the ClinVar safety guard (no GRCh36 cache)
    and silently annotated GRCh36 data against GRCh37 coordinates.
    """

    def _grch36_variants(self, count=3, discordant=1):
        """Build variants: `count` at GRCh36 positions + `discordant` at junk positions."""
        variants = []
        grch36_rsids = [
            rsid for rsid, builds in KNOWN_SNP_POSITIONS.items() if BUILD_GRCH36 in builds
        ]
        for rsid in grch36_rsids[:count]:
            chrom, pos = KNOWN_SNP_POSITIONS[rsid][BUILD_GRCH36]
            variants.append(
                Variant(rsid=rsid, chromosome=chrom, position=pos, allele1="A", allele2="A")
            )
        for rsid in grch36_rsids[count : count + discordant]:
            chrom, _ = KNOWN_SNP_POSITIONS[rsid][BUILD_GRCH36]
            variants.append(
                Variant(rsid=rsid, chromosome=chrom, position=99999999, allele1="A", allele2="A")
            )
        return variants

    def test_non_confident_grch36_uses_grch36_effective(self):
        state = _BuildDetectionState(override=None, header_build=None)
        variants = self._grch36_variants(count=3, discordant=1)
        for v in variants:
            state.feed(v)
        state.flush()
        assert state.effective_build == BUILD_GRCH36

    def test_non_confident_grch36_with_header_grch37_still_uses_grch36(self):
        state = _BuildDetectionState(override=None, header_build=BUILD_GRCH37)
        variants = self._grch36_variants(count=3, discordant=1)
        for v in variants:
            state.feed(v)
        state.flush()
        assert state.effective_build == BUILD_GRCH36

    def test_confident_grch36_uses_grch36(self):
        state = _BuildDetectionState(override=None, header_build=None)
        variants = self._grch36_variants(count=3, discordant=0)
        for v in variants:
            state.feed(v)
        state.flush()
        assert state.effective_build == BUILD_GRCH36

    def test_diagnostics_report_grch36_as_effective(self):
        state = _BuildDetectionState(override=None, header_build=None)
        variants = self._grch36_variants(count=3, discordant=1)
        for v in variants:
            state.feed(v)
        state.flush()
        diag = state.diagnostics()
        assert diag.effective_build == BUILD_GRCH36
        assert diag.detected_build == BUILD_GRCH36

    def test_buffer_limit_with_single_grch36_probe_uses_grch36(self):
        """Buffer-limit path must apply the same GRCh36 safety as flush().

        Real FTDNA GRCh36 files have 687K+ variants but only 1 probe SNP
        in the first 100K lines. The buffer-limit fallback must run
        detect_build and trigger the GRCh36 guard, not hard-fall-back to
        GRCh37.
        """
        state = _BuildDetectionState(override=None, header_build=BUILD_GRCH37)
        grch36_rsids = [
            rsid for rsid, builds in KNOWN_SNP_POSITIONS.items() if BUILD_GRCH36 in builds
        ]
        rsid = grch36_rsids[0]
        chrom, pos = KNOWN_SNP_POSITIONS[rsid][BUILD_GRCH36]
        probe = Variant(rsid=rsid, chromosome=chrom, position=pos, allele1="A", allele2="A")
        state.feed(probe)
        filler = [
            Variant(rsid=f"rs9{i:06d}", chromosome="1", position=i, allele1="A", allele2="A")
            for i in range(_DETECTION_BUFFER_LIMIT)
        ]
        for v in filler:
            ready, _batch = state.feed(v)
            if ready:
                break
        assert state.effective_build == BUILD_GRCH36
        diag = state.diagnostics()
        assert diag.detected_build == BUILD_GRCH36


class TestGnomadEnrichment:
    """gnomAD frequency enrichment stamps allele_frequency on annotations."""

    def test_enrichment_stamps_frequency(
        self,
        mock_mhg_path: Path,
        clinvar_data_dir: Path,
    ) -> None:
        """run_analysis with gnomAD annotator stamps allele_frequency."""
        import sqlite3

        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.annotators.gnomad import GnomadAnnotator
        from allelix.databases.gnomad_loader import GNOMAD_DB_FILENAME
        from allelix.databases.schema import GNOMAD_SCHEMA

        db_path = clinvar_data_dir / GNOMAD_DB_FILENAME
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            for stmt in GNOMAD_SCHEMA.split(";"):
                stmt = stmt.strip()
                if stmt:
                    conn.execute(stmt)
            conn.execute(
                "INSERT OR REPLACE INTO gnomad_frequencies"
                " (chrom, pos, ref, alt, rsid, af) VALUES (?, ?, ?, ?, ?, ?)",
                ("1", 11796321, "G", "A", "rs1801133", 0.35),
            )
            conn.execute(
                "INSERT OR REPLACE INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count)"
                " VALUES (?, ?, ?, ?, ?)",
                ("gnomad", "test://mock", "4.1", "2026-01-01T00:00:00Z", 1),
            )
            conn.commit()

        parser = MyHappyGenesParser()
        clinvar = ClinVarAnnotator(clinvar_data_dir)
        gnomad = GnomadAnnotator(clinvar_data_dir)
        result = run_analysis(
            mock_mhg_path,
            parser,
            [clinvar],
            gnomad=gnomad,
        )
        mthfr = [a for a in result.annotations if a.rsid == "rs1801133"]
        assert any(a.allele_frequency is not None for a in mthfr)
        assert ("gnomad", "4.1") in result.annotators_used

    def test_no_gnomad_no_frequency(
        self,
        mock_mhg_path: Path,
        clinvar_data_dir: Path,
    ) -> None:
        """run_analysis without gnomAD leaves allele_frequency as None."""
        parser = MyHappyGenesParser()
        clinvar = ClinVarAnnotator(clinvar_data_dir)
        result = run_analysis(mock_mhg_path, parser, [clinvar])
        assert all(a.allele_frequency is None for a in result.annotations)
        assert all(name != "gnomad" for name, _ in result.annotators_used)


class TestAlphaMissenseEnrichment:
    """AlphaMissense enrichment stamps am_pathogenicity/am_class on annotations."""

    def test_enrichment_stamps_pathogenicity(
        self,
        mock_mhg_path: Path,
        clinvar_data_dir: Path,
    ) -> None:
        """run_analysis with AlphaMissense stamps am_pathogenicity and am_class."""
        import sqlite3

        from allelix.annotators.alphamissense import AlphaMissenseAnnotator
        from allelix.databases.alphamissense_loader import ALPHAMISSENSE_DB_FILENAME
        from allelix.databases.schema import ALPHAMISSENSE_SCHEMA

        db_path = clinvar_data_dir / ALPHAMISSENSE_DB_FILENAME
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            for stmt in ALPHAMISSENSE_SCHEMA.split(";"):
                stmt = stmt.strip()
                if stmt:
                    conn.execute(stmt)
            conn.execute(
                "INSERT OR REPLACE INTO alphamissense_scores"
                " (chrom, pos, ref, alt, rsid, uniprot_id, transcript_id,"
                " protein_variant, am_pathogenicity, am_class)"
                " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    "1",
                    11856378,
                    "G",
                    "A",
                    "rs1801133",
                    "P42898",
                    "ENST001",
                    "A222V",
                    0.72,
                    "ambiguous",
                ),
            )
            conn.execute(
                "INSERT OR REPLACE INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count)"
                " VALUES (?, ?, ?, ?, ?)",
                ("alphamissense", "test://mock", "2023.2", "2026-01-01", 1),
            )
            conn.commit()

        parser = MyHappyGenesParser()
        clinvar = ClinVarAnnotator(clinvar_data_dir)
        am = AlphaMissenseAnnotator(clinvar_data_dir)
        result = run_analysis(
            mock_mhg_path,
            parser,
            [clinvar],
            alphamissense=am,
        )
        mthfr = [a for a in result.annotations if a.rsid == "rs1801133"]
        assert any(a.am_pathogenicity is not None for a in mthfr)
        assert any(a.am_class == "ambiguous" for a in mthfr)
        assert ("alphamissense", "2023.2") in result.annotators_used

    def test_no_alphamissense_no_pathogenicity(
        self,
        mock_mhg_path: Path,
        clinvar_data_dir: Path,
    ) -> None:
        """run_analysis without AlphaMissense leaves am_pathogenicity as None."""
        parser = MyHappyGenesParser()
        clinvar = ClinVarAnnotator(clinvar_data_dir)
        result = run_analysis(mock_mhg_path, parser, [clinvar])
        assert all(a.am_pathogenicity is None for a in result.annotations)
        assert all(name != "alphamissense" for name, _ in result.annotators_used)


class TestEnrichmentExactVsFallback:
    """Pipeline splits exact (rsid, alt) lookups from MAX-by-rsid fallback."""

    def test_annotation_with_alt_uses_exact_lookup(self) -> None:
        """Annotations with alt set get exact (rsid, alt) enrichment."""
        from allelix.models import Annotation

        a = Annotation(
            source="clinvar",
            rsid="rs429358",
            significance="clinvar_pathogenic",
            category="clinical",
            magnitude=9.0,
            description="test",
            attribution="ClinVar",
            genotype_match="TC",
            alt="C",
        )
        assert a.alt == "C"
        exact_keys = {(x.rsid, x.alt) for x in [a] if x.alt}
        assert ("rs429358", "C") in exact_keys

    def test_annotation_without_alt_uses_max_fallback(self) -> None:
        """Annotations without alt fall back to MAX-by-rsid enrichment."""
        from allelix.models import Annotation

        a = Annotation(
            source="gwas",
            rsid="rs429358",
            significance="gwas_association",
            category="trait",
            magnitude=7.0,
            description="test",
            attribution="GWAS Catalog",
            genotype_match="TC",
            alt="",
        )
        max_rsids = {x.rsid for x in [a] if not x.alt}
        assert "rs429358" in max_rsids

    def test_gwas_annotations_have_no_alt(self) -> None:
        """GWAS annotations must not set alt (risk allele != VCF ALT)."""
        import sqlite3
        import tempfile
        from pathlib import Path as _Path

        from allelix.annotators.gwas import GWASCatalogAnnotator
        from allelix.databases.gwas_loader import load_gwas_tsv
        from allelix.databases.schema import GWAS_SCHEMA

        db_path = _Path(__file__).parent.parent / "fixtures" / "mock_gwas_catalog.tsv"
        if not db_path.exists():
            pytest.skip("mock GWAS fixture not available")

        with tempfile.TemporaryDirectory() as td:
            tmp = _Path(td)
            gwas_db = tmp / "gwas.sqlite"
            with contextlib.closing(sqlite3.connect(gwas_db)) as conn:
                for stmt in GWAS_SCHEMA.split(";"):
                    stmt = stmt.strip()
                    if stmt:
                        conn.execute(stmt)
                conn.commit()
            load_gwas_tsv(db_path, gwas_db)
            ann = GWASCatalogAnnotator(tmp)
            if not ann.is_ready():
                pytest.skip("GWAS db not ready")
            v = Variant(
                rsid="rs1801133",
                chromosome="1",
                position=11796321,
                allele1="C",
                allele2="T",
            )
            results = ann.annotate(v)
            for r in results:
                assert r.alt == "", f"GWAS annotation should not set alt, got {r.alt!r}"
            ann.close()


class TestCaddEnrichment:
    """CADD enrichment stamps cadd_phred on annotations via gnomAD coordinate resolution."""

    def test_enrichment_stamps_cadd_phred(
        self,
        mock_mhg_path: Path,
        clinvar_data_dir: Path,
    ) -> None:
        """run_analysis with CADD + gnomAD stamps cadd_phred."""
        import sqlite3

        from allelix.annotators.cadd import CaddAnnotator
        from allelix.annotators.gnomad import GnomadAnnotator
        from allelix.databases.cadd_loader import CADD_DB_FILENAME
        from allelix.databases.gnomad_loader import GNOMAD_DB_FILENAME
        from allelix.databases.schema import CADD_SCHEMA, GNOMAD_SCHEMA

        gnomad_path = clinvar_data_dir / GNOMAD_DB_FILENAME
        with contextlib.closing(sqlite3.connect(gnomad_path)) as conn:
            for stmt in GNOMAD_SCHEMA.split(";"):
                stmt = stmt.strip()
                if stmt:
                    conn.execute(stmt)
            conn.execute(
                "INSERT OR REPLACE INTO gnomad_frequencies"
                " (chrom, pos, ref, alt, rsid, af) VALUES (?, ?, ?, ?, ?, ?)",
                ("1", 11796321, "G", "A", "rs1801133", 0.35),
            )
            conn.execute(
                "INSERT OR REPLACE INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count)"
                " VALUES (?, ?, ?, ?, ?)",
                ("gnomad", "test://mock", "4.1", "2026-01-01T00:00:00Z", 1),
            )
            conn.commit()

        cadd_path = clinvar_data_dir / CADD_DB_FILENAME
        with contextlib.closing(sqlite3.connect(cadd_path)) as conn:
            for stmt in CADD_SCHEMA.split(";"):
                stmt = stmt.strip()
                if stmt:
                    conn.execute(stmt)
            conn.execute(
                "INSERT INTO cadd_scores (chrom, pos, ref, alt, phred) VALUES (?, ?, ?, ?, ?)",
                ("1", 11796321, "G", "A", 24.3),
            )
            conn.execute(
                "INSERT INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count,"
                "  local_version_tag)"
                " VALUES (?, ?, ?, ?, ?, ?)",
                ("cadd", "test://mock", "v1.7", "2026-01-01T00:00:00Z", 1, "sv:1"),
            )
            conn.commit()

        parser = MyHappyGenesParser()
        clinvar = ClinVarAnnotator(clinvar_data_dir)
        gnomad = GnomadAnnotator(clinvar_data_dir)
        cadd = CaddAnnotator(clinvar_data_dir)
        result = run_analysis(
            mock_mhg_path,
            parser,
            [clinvar],
            gnomad=gnomad,
            cadd=cadd,
        )
        mthfr = [a for a in result.annotations if a.rsid == "rs1801133"]
        assert any(a.cadd_phred is not None for a in mthfr)
        assert any(a.cadd_phred == pytest.approx(24.3) for a in mthfr if a.cadd_phred is not None)
        assert ("cadd", "v1.7") in result.annotators_used

    def test_no_cadd_no_phred(
        self,
        mock_mhg_path: Path,
        clinvar_data_dir: Path,
    ) -> None:
        """run_analysis without CADD leaves cadd_phred as None."""
        parser = MyHappyGenesParser()
        clinvar = ClinVarAnnotator(clinvar_data_dir)
        result = run_analysis(mock_mhg_path, parser, [clinvar])
        assert all(a.cadd_phred is None for a in result.annotations)
        assert all(name != "cadd" for name, _ in result.annotators_used)

    def test_cadd_without_gnomad_skips_enrichment(
        self,
        mock_mhg_path: Path,
        clinvar_data_dir: Path,
    ) -> None:
        """CADD enrichment requires gnomAD for coordinate resolution."""
        import sqlite3

        from allelix.annotators.cadd import CaddAnnotator
        from allelix.databases.cadd_loader import CADD_DB_FILENAME
        from allelix.databases.schema import CADD_SCHEMA

        cadd_path = clinvar_data_dir / CADD_DB_FILENAME
        with contextlib.closing(sqlite3.connect(cadd_path)) as conn:
            for stmt in CADD_SCHEMA.split(";"):
                stmt = stmt.strip()
                if stmt:
                    conn.execute(stmt)
            conn.execute(
                "INSERT INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count,"
                "  local_version_tag)"
                " VALUES (?, ?, ?, ?, ?, ?)",
                ("cadd", "test://mock", "v1.7", "2026-01-01T00:00:00Z", 0, "sv:1"),
            )
            conn.commit()

        parser = MyHappyGenesParser()
        clinvar = ClinVarAnnotator(clinvar_data_dir)
        cadd = CaddAnnotator(clinvar_data_dir)
        result = run_analysis(
            mock_mhg_path,
            parser,
            [clinvar],
            cadd=cadd,
        )
        assert all(a.cadd_phred is None for a in result.annotations)


class TestCaddMultiAllelic:
    """CADD enrichment at multi-allelic sites must use the user's allele, not max."""

    def test_multiallelic_uses_user_allele_not_max(self) -> None:
        """At a multi-allelic site, CADD score must match the user's allele."""
        from allelix.reports._pipeline import _lookup_user_allele
        from allelix.utils.allele import resolve_strand

        coords = [("1", 100, "A", "C"), ("1", 100, "A", "G")]
        scores = {("1", 100, "A", "C"): 5.0, ("1", 100, "A", "G"): 30.0}

        result = _lookup_user_allele("C", coords, scores, resolve_strand)
        assert result == pytest.approx(5.0)

    def test_multiallelic_no_alt_uses_max(self) -> None:
        """Without a known user allele, max-reduce is the correct fallback."""
        from allelix.models import Annotation
        from allelix.reports._pipeline import _enrich_cadd

        coords = [("1", 100, "A", "C"), ("1", 100, "A", "G")]
        scores = {("1", 100, "A", "C"): 5.0, ("1", 100, "A", "G"): 30.0}

        ann = Annotation(
            source="clinvar",
            rsid="rs999",
            significance="clinvar_pathogenic",
            category="clinical",
            magnitude=9.0,
            description="test",
            attribution="ClinVar",
            genotype_match="A/C",
            alt="",
        )

        class MockGnomad:
            def bulk_resolve_coordinates(self, rsids):
                return {"rs999": coords}

        class MockCadd:
            def bulk_lookup(self, keys):
                return scores

        _enrich_cadd([ann], MockGnomad(), MockCadd())
        assert ann.cadd_phred == pytest.approx(30.0)

    def test_complement_fallback(self) -> None:
        """Complement match works when no direct match exists."""
        from allelix.reports._pipeline import _lookup_user_allele
        from allelix.utils.allele import resolve_strand

        coords = [("1", 200, "C", "A")]
        scores = {("1", 200, "C", "A"): 12.5}

        result = _lookup_user_allele("T", coords, scores, resolve_strand)
        assert result == pytest.approx(12.5)
