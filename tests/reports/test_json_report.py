# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Tests for the JSON report renderer."""

from __future__ import annotations

import json
from pathlib import Path

from allelix.models import Annotation
from allelix.reports._pipeline import AnalysisResult
from allelix.reports.json_report import REGULATORY_NOTICE, SCHEMA_VERSION, render_json


def _result(annotations: list[Annotation]) -> AnalysisResult:
    return AnalysisResult(
        file_path=Path("genotype.txt"),
        parser_name="myhappygenes",
        parser_display_name="MyHappyGenes (Tempus)",
        sample_id="MHG_X",
        build="GRCh37",
        total_variants=10,
        skipped_count=0,
        annotators_used=[("clinvar", "20260101"), ("pharmgkb", "2026-05-11")],
        annotations=annotations,
    )


def _ann(**overrides) -> Annotation:
    defaults = {
        "source": "clinvar",
        "rsid": "rs1801133",
        "significance": "clinvar_pathogenic",
        "category": "clinical",
        "magnitude": 9.0,
        "description": "x",
        "attribution": "ClinVar",
        "genotype_match": "A",
        "gene": "MTHFR",
        "condition": "MTHFR deficiency",
        "references": ["clinvar:allele/100001"],
    }
    defaults.update(overrides)
    return Annotation(**defaults)


class TestRenderJson:
    def test_writes_valid_json(self, tmp_path: Path):
        out = tmp_path / "report.json"
        count = render_json(_result([_ann()]), output_path=out)
        assert count == 1
        payload = json.loads(out.read_text())
        assert payload["schema_version"] == SCHEMA_VERSION
        assert payload["input"]["sample_id"] == "MHG_X"
        assert payload["input"]["total_variants"] == 10
        assert payload["regulatory_notice"] == REGULATORY_NOTICE

    def test_annotators_listed_with_versions(self, tmp_path: Path):
        out = tmp_path / "r.json"
        render_json(_result([_ann()]), output_path=out)
        payload = json.loads(out.read_text())
        assert {"name": "clinvar", "version": "20260101"} in payload["annotators"]
        assert {"name": "pharmgkb", "version": "2026-05-11"} in payload["annotators"]

    def test_annotation_round_trip(self, tmp_path: Path):
        out = tmp_path / "r.json"
        render_json(_result([_ann(rsid="rs99")]), output_path=out)
        payload = json.loads(out.read_text())
        ann = payload["annotations"][0]
        assert ann["rsid"] == "rs99"
        assert ann["attribution"] == "ClinVar"  # ADR-0003: every row attributed
        assert ann["significance"].startswith("clinvar_")

    def test_filters_recorded_in_payload(self, tmp_path: Path):
        out = tmp_path / "r.json"
        render_json(
            _result([_ann()]),
            output_path=out,
            min_magnitude=5.0,
            category="clinical",
            genes={"MTHFR"},
        )
        payload = json.loads(out.read_text())
        assert payload["filters"]["min_magnitude"] == 5.0
        assert payload["filters"]["category"] == "clinical"
        assert payload["filters"]["genes"] == ["MTHFR"]

    def test_min_magnitude_filters_annotations(self, tmp_path: Path):
        anns = [_ann(rsid="lo", magnitude=2.0), _ann(rsid="hi", magnitude=8.0)]
        out = tmp_path / "r.json"
        count = render_json(_result(anns), output_path=out, min_magnitude=5.0)
        assert count == 1
        payload = json.loads(out.read_text())
        assert [a["rsid"] for a in payload["annotations"]] == ["hi"]

    def test_empty_annotations(self, tmp_path: Path):
        out = tmp_path / "r.json"
        count = render_json(_result([]), output_path=out)
        assert count == 0
        payload = json.loads(out.read_text())
        assert payload["annotations"] == []
