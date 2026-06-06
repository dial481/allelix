# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Tests for terminal report rendering."""

from __future__ import annotations

import io
from pathlib import Path

from rich.console import Console

from allelix.models import Annotation
from allelix.reports._pipeline import AnalysisResult
from allelix.reports.terminal import render_terminal, render_terminal_diff


def _ann(**overrides) -> Annotation:
    defaults = {
        "source": "clinvar",
        "rsid": "rs1",
        "significance": "clinvar_pathogenic",
        "category": "clinical",
        "magnitude": 5.0,
        "description": "test",
        "attribution": "ClinVar",
        "genotype_match": "A",
        "gene": "GENE1",
        "condition": "Some condition",
    }
    defaults.update(overrides)
    return Annotation(**defaults)


def _result(annotations: list[Annotation]) -> AnalysisResult:
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


def _render(annotations, **kwargs) -> tuple[str, int]:
    buf = io.StringIO()
    console = Console(file=buf, force_terminal=False, width=200)
    count = render_terminal(_result(annotations), console=console, **kwargs)
    return buf.getvalue(), count


class TestRenderTerminal:
    def test_empty_list_renders_message(self):
        out, count = _render([])
        assert count == 0
        assert "No annotations" in out

    def test_attribution_column_present(self):
        out, _ = _render([_ann()])
        assert "ClinVar" in out
        assert "Source" in out

    def test_sorts_by_magnitude_descending(self):
        annotations = [
            _ann(rsid="rs_low", magnitude=2.0),
            _ann(rsid="rs_high", magnitude=9.0),
            _ann(rsid="rs_mid", magnitude=5.0),
        ]
        out, _ = _render(annotations)
        assert out.index("rs_high") < out.index("rs_mid") < out.index("rs_low")

    def test_min_magnitude_filter(self):
        annotations = [
            _ann(rsid="rs_skip", magnitude=2.0),
            _ann(rsid="rs_keep", magnitude=8.0),
        ]
        out, count = _render(annotations, min_magnitude=5.0)
        assert count == 1
        assert "rs_keep" in out
        assert "rs_skip" not in out

    def test_category_filter(self):
        annotations = [
            _ann(rsid="rs_clinical", category="clinical"),
            _ann(rsid="rs_pharma", category="pharma"),
        ]
        out, _ = _render(annotations, category="clinical")
        assert "rs_clinical" in out
        assert "rs_pharma" not in out

    def test_genes_filter(self):
        annotations = [
            _ann(rsid="m", gene="MTHFR"),
            _ann(rsid="b", gene="BRCA1"),
        ]
        out, count = _render(annotations, genes={"MTHFR"})
        assert count == 1
        assert "m" in out
        assert "BRCA1" not in out


def _ann_dict(**overrides) -> dict:
    defaults = {
        "source": "clinvar",
        "rsid": "rs1801133",
        "significance": "clinvar_pathogenic",
        "category": "clinical",
        "magnitude": 9.0,
        "description": "clinvar: test",
        "attribution": "ClinVar",
        "genotype_match": "AG",
        "references": [],
        "condition": "MTHFR deficiency",
        "gene": "MTHFR",
    }
    defaults.update(overrides)
    return defaults


class TestRenderTerminalDiff:
    def test_render_terminal_diff_new_only(self, capsys):
        """Verifies the New Annotations table renders with all columns."""
        from allelix.reports.diff import DiffResult

        diff = DiffResult(
            new=[_ann(rsid="rs1801133", gene="MTHFR", magnitude=9.0)],
            previous_generated_at="2026-05-01T00:00:00",
        )
        total = render_terminal_diff(diff, Console(force_terminal=True, width=200))
        out = capsys.readouterr().out
        assert total == 1
        assert "New Annotations (1)" in out
        assert "rs1801133" in out

    def test_render_terminal_diff_changed_only(self, capsys):
        """Old Sig / New Sig / Old Mag / New Mag columns all present."""
        from allelix.reports.diff import ChangedAnnotation, DiffResult

        diff = DiffResult(
            changed=[
                ChangedAnnotation(
                    current=_ann(magnitude=7.0),
                    previous_significance="old_sig",
                    previous_magnitude=9.0,
                )
            ],
            previous_generated_at="2026-05-01T00:00:00",
        )
        render_terminal_diff(diff, Console(force_terminal=True, width=200))
        out = capsys.readouterr().out
        assert "Changed Annotations (1)" in out
        assert "old_sig" in out
        assert "9.0" in out and "7.0" in out

    def test_render_terminal_diff_removed_only(self, capsys):
        from allelix.reports.diff import DiffResult

        diff = DiffResult(
            removed=[_ann_dict()],
            previous_generated_at="2026-05-01T00:00:00",
        )
        render_terminal_diff(diff, Console(force_terminal=True, width=200))
        out = capsys.readouterr().out
        assert "Removed Annotations (1)" in out
        assert "rs1801133" in out

    def test_render_terminal_diff_no_changes(self, capsys):
        from allelix.reports.diff import DiffResult

        diff = DiffResult(previous_generated_at="2026-05-01T00:00:00")
        total = render_terminal_diff(diff, Console(force_terminal=True, width=200))
        assert total == 0
        assert "No changes since previous report." in capsys.readouterr().out
