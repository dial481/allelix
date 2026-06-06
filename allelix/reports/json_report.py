# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""JSON report renderer.

Output schema (versioned via `schema_version`):

    {
      "schema_version": "1",
      "allelix_version": "0.4.0",
      "generated_at": "2026-05-11T12:34:56+00:00",
      "regulatory_notice": "...",
      "input": {
        "file": "genotype.txt",
        "format": "myhappygenes",
        "sample_id": "MHG000001",
        "build": "GRCh37",
        "total_variants": 2015,
        "skipped_lines": 0
      },
      "annotators": [
        {"name": "clinvar", "version": "20260101"}
      ],
      "filters": {
        "min_magnitude": 5.0,
        "category": null,
        "genes": null
      },
      "annotations": [ ... ]
    }

Every annotation is source-attributed (ADR-0003); the renderer never adds
or strips that field.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from datetime import UTC, datetime
from typing import TYPE_CHECKING

from allelix import __version__
from allelix.reports import REGULATORY_NOTICE, atomic_write_text
from allelix.reports._pipeline import rollup_gwas_duplicates

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path

    from allelix.reports._pipeline import AnalysisResult
    from allelix.reports.diff import DiffResult


SCHEMA_VERSION = "1"

__all__ = ["REGULATORY_NOTICE", "SCHEMA_VERSION", "render_json"]


def render_json(
    result: AnalysisResult,
    *,
    output_path: Path,
    min_magnitude: float = 0.0,
    category: str | None = None,
    genes: Iterable[str] | None = None,
    source_min_magnitudes: dict[str, float] | None = None,
    diff: DiffResult | None = None,
) -> int:
    """Write a JSON report to `output_path`. Returns the number of annotations included."""
    filtered = result.filter(
        min_magnitude=min_magnitude,
        category=category,
        genes=genes,
        source_min_magnitudes=source_min_magnitudes,
    )
    filtered = rollup_gwas_duplicates(filtered)
    payload: dict = {
        "schema_version": SCHEMA_VERSION,
        "allelix_version": __version__,
        "generated_at": datetime.now(UTC).isoformat(),
        "regulatory_notice": REGULATORY_NOTICE,
        "input": {
            "file": result.file_path.name,
            "format": result.parser_name,
            "sample_id": result.sample_id,
            "build": result.build,
            "total_variants": result.total_variants,
            "skipped_lines": result.skipped_count,
        },
        "annotators": [
            {"name": name, "version": version} for name, version in result.annotators_used
        ],
        "filters": {
            "min_magnitude": min_magnitude,
            "category": category,
            "genes": sorted(genes) if genes else None,
        },
        "annotations": [asdict(a) for a in filtered],
    }

    if diff is not None:
        from allelix.reports.diff import diff_annotation_to_dict, summarize_diff

        payload["diff"] = {
            "previous_report": diff.previous_generated_at,
            "summary": summarize_diff(diff),
            "new": [asdict(a) for a in diff.new],
            "changed": [diff_annotation_to_dict(c) for c in diff.changed],
            "removed": diff.removed,
        }

    atomic_write_text(output_path, json.dumps(payload, indent=2, sort_keys=False) + "\n")
    return len(filtered)
