# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Tests for the unified analysis pipeline."""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.parsers.myhappygenes import MyHappyGenesParser
from allelix.reports._pipeline import AnalysisResult, run_analysis

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
