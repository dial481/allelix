# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Tests for the HTML report renderer."""

from __future__ import annotations

from pathlib import Path

from allelix.models import Annotation
from allelix.reports._pipeline import AnalysisResult, BuildDiagnostics
from allelix.reports.html import render_html


def _result(annotations: list[Annotation]) -> AnalysisResult:
    return AnalysisResult(
        file_path=Path("genotype.txt"),
        parser_name="myhappygenes",
        parser_display_name="MyHappyGenes (Tempus)",
        sample_id="MHG_HTML",
        build="GRCh37",
        total_variants=10,
        skipped_count=0,
        annotators_used=[("clinvar", "20260101")],
        annotations=annotations,
    )


def _ann(**overrides) -> Annotation:
    defaults = {
        "source": "clinvar",
        "rsid": "rs1801133",
        "significance": "clinvar_pathogenic",
        "category": "clinical",
        "magnitude": 9.0,
        "description": "ClinVar classifies this allele as Pathogenic",
        "attribution": "ClinVar",
        "genotype_match": "A",
        "gene": "MTHFR",
        "condition": "MTHFR deficiency",
    }
    defaults.update(overrides)
    return Annotation(**defaults)


class TestRenderHtml:
    def test_writes_self_contained_file(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert body.startswith("<!DOCTYPE html>")
        # Inline CSS — no external stylesheet links
        assert "<style>" in body
        assert '<link rel="stylesheet"' not in body
        assert "<script src=" not in body

    def test_disclaimer_present(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert "Informational only" in body
        assert "not medical advice" in body  # from REGULATORY_NOTICE

    def test_attribution_visible_in_table(self, tmp_path: Path):
        """ADR-0003: every row's source attribution must render."""
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert "ClinVar" in body
        assert "rs1801133" in body
        assert "MTHFR" in body

    def test_html_escapes_user_supplied_strings(self, tmp_path: Path):
        """A condition string with HTML must be escaped, not rendered as markup."""
        evil = _ann(condition='<script>alert("xss")</script>')
        out = tmp_path / "report.html"
        render_html(_result([evil]), output_path=out)
        body = out.read_text()
        assert "<script>alert" not in body
        assert "&lt;script&gt;" in body

    def test_empty_annotations(self, tmp_path: Path):
        out = tmp_path / "report.html"
        count = render_html(_result([]), output_path=out)
        assert count == 0
        body = out.read_text()
        assert "No annotations" in body

    def test_min_magnitude_filter(self, tmp_path: Path):
        anns = [_ann(rsid="lo", magnitude=2.0), _ann(rsid="hi", magnitude=8.0)]
        out = tmp_path / "report.html"
        count = render_html(_result(anns), output_path=out, min_magnitude=5.0)
        assert count == 1
        body = out.read_text()
        assert "rs_hi" not in body  # we use "hi"
        assert '<td class="rsid">hi</td>' in body
        assert '<td class="rsid">lo</td>' not in body

    def test_genes_filter(self, tmp_path: Path):
        anns = [_ann(rsid="m", gene="MTHFR"), _ann(rsid="b", gene="BRCA1")]
        out = tmp_path / "report.html"
        render_html(_result(anns), output_path=out, genes={"MTHFR"})
        body = out.read_text()
        assert ">m<" in body
        assert ">b<" not in body


class TestEducationSection:
    def test_education_section_present(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert "Reading This Report" in body

    def test_education_pseudogene_warning(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert "cross-hybridize" in body

    def test_education_carrier_vs_affected(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert "Carrier vs. affected" in body

    def test_education_confirmatory_testing(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert "Confirmatory testing" in body


class TestBuildMismatchBanner:
    def test_no_banner_without_mismatch(self, tmp_path: Path):
        r = _result([_ann()])
        r.build_diagnostics = BuildDiagnostics(
            header_build="GRCh37",
            detected_build="GRCh37",
            effective_build="GRCh37",
            override=False,
            matched_count=5,
            inspected_count=5,
        )
        out = tmp_path / "report.html"
        render_html(r, output_path=out)
        body = out.read_text()
        assert "Build mismatch" not in body

    def test_banner_on_mismatch(self, tmp_path: Path):
        r = _result([_ann()])
        r.build_diagnostics = BuildDiagnostics(
            header_build="GRCh37",
            detected_build="GRCh38",
            effective_build="GRCh38",
            override=False,
            matched_count=5,
            inspected_count=5,
        )
        out = tmp_path / "report.html"
        render_html(r, output_path=out)
        body = out.read_text()
        assert "Build mismatch" in body
        assert "GRCh37" in body
        assert "GRCh38" in body
        assert "notice-warn" in body

    def test_no_banner_with_override(self, tmp_path: Path):
        r = _result([_ann()])
        r.build_diagnostics = BuildDiagnostics(
            header_build="GRCh37",
            detected_build="GRCh38",
            effective_build="GRCh38",
            override=True,
            matched_count=5,
            inspected_count=5,
        )
        out = tmp_path / "report.html"
        render_html(r, output_path=out)
        body = out.read_text()
        assert "Build mismatch" not in body

    def test_no_banner_without_diagnostics(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert "Build mismatch" not in body
