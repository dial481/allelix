# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Terminal report rendering for `allelix analyze`."""

from __future__ import annotations

from typing import TYPE_CHECKING

from rich.table import Table

from allelix.reports._pipeline import rollup_gwas_duplicates

if TYPE_CHECKING:
    from collections.abc import Iterable

    from rich.console import Console

    from allelix.models import Annotation
    from allelix.reports._pipeline import AnalysisResult
    from allelix.reports.diff import DiffResult


def render_terminal(
    result: AnalysisResult,
    console: Console,
    *,
    min_magnitude: float = 0.0,
    category: str | None = None,
    genes: Iterable[str] | None = None,
    source_min_magnitudes: dict[str, float] | None = None,
) -> int:
    """Render an AnalysisResult as a Rich table. Returns annotation count.

    Per ADR-0003 (regulatory posture), every row shows the source attribution
    in its own column — no rendered claim is unattributed.
    """
    filtered = result.filter(
        min_magnitude=min_magnitude,
        category=category,
        genes=genes,
        source_min_magnitudes=source_min_magnitudes,
    )
    filtered = rollup_gwas_duplicates(filtered)
    _print_table(filtered, console)
    return len(filtered)


def render_terminal_diff(
    diff: DiffResult,
    console: Console,
) -> int:
    """Render a diff summary and tables for new/changed/removed annotations."""
    from allelix.reports.diff import summarize_diff

    summary = summarize_diff(diff)
    if not diff.has_changes:
        console.print(f"[green]{summary}[/green]")
        return 0

    console.print(f"[bold]{summary}[/bold]")
    total = 0

    if diff.new:
        table = Table(title=f"New Annotations ({len(diff.new)})")
        table.add_column("rsID", style="cyan", no_wrap=True)
        table.add_column("Gene", style="magenta", no_wrap=True)
        table.add_column("Source", style="blue", no_wrap=True)
        table.add_column("Significance", style="yellow")
        table.add_column("Review Status", style="dim")
        table.add_column("Magnitude", justify="right")
        table.add_column("Genotype", no_wrap=True)
        table.add_column("Condition", overflow="fold")
        for a in diff.new:
            table.add_row(
                a.rsid,
                a.gene or "—",
                a.attribution,
                a.significance,
                a.review_status or "—",
                f"{a.magnitude:.1f}",
                a.genotype_match,
                a.condition or "—",
            )
        console.print(table)
        total += len(diff.new)

    if diff.changed:
        table = Table(title=f"Changed Annotations ({len(diff.changed)})")
        table.add_column("rsID", style="cyan", no_wrap=True)
        table.add_column("Gene", style="magenta", no_wrap=True)
        table.add_column("Source", style="blue", no_wrap=True)
        table.add_column("Old Sig", style="dim")
        table.add_column("New Sig", style="yellow")
        table.add_column("Review Status", style="dim")
        table.add_column("Old Mag", justify="right", style="dim")
        table.add_column("New Mag", justify="right")
        table.add_column("Condition", overflow="fold")
        for c in diff.changed:
            table.add_row(
                c.current.rsid,
                c.current.gene or "—",
                c.current.attribution,
                c.previous_significance,
                c.current.significance,
                c.current.review_status or "—",
                f"{c.previous_magnitude:.1f}",
                f"{c.current.magnitude:.1f}",
                c.current.condition or "—",
            )
        console.print(table)
        total += len(diff.changed)

    if diff.removed:
        table = Table(title=f"Removed Annotations ({len(diff.removed)})")
        table.add_column("rsID", style="dim cyan", no_wrap=True)
        table.add_column("Gene", style="dim magenta", no_wrap=True)
        table.add_column("Source", style="dim blue", no_wrap=True)
        table.add_column("Significance", style="dim")
        table.add_column("Review Status", style="dim")
        table.add_column("Magnitude", justify="right", style="dim")
        table.add_column("Condition", overflow="fold", style="dim")
        for d in diff.removed:
            table.add_row(
                d.get("rsid", ""),
                d.get("gene", "") or "—",
                d.get("attribution", ""),
                d.get("significance", ""),
                d.get("review_status", "") or "—",
                f"{d.get('magnitude', 0.0):.1f}",
                d.get("condition", "") or "—",
            )
        console.print(table)
        total += len(diff.removed)

    return total


def _print_table(filtered: list[Annotation], console: Console) -> None:
    if not filtered:
        console.print("[yellow]No annotations matched the current filters.[/yellow]")
        return

    table = Table(title=f"Annotations ({len(filtered)})")
    table.add_column("rsID", style="cyan", no_wrap=True)
    table.add_column("Gene", style="magenta", no_wrap=True)
    table.add_column("Source", style="blue", no_wrap=True)
    table.add_column("Significance", style="yellow")
    table.add_column("Review Status", style="dim")
    table.add_column("Magnitude", justify="right")
    table.add_column("Genotype", no_wrap=True)
    table.add_column("Condition", overflow="fold")

    for a in filtered:
        table.add_row(
            a.rsid,
            a.gene or "—",
            a.attribution,
            a.significance,
            a.review_status or "—",
            f"{a.magnitude:.1f}",
            a.genotype_match,
            a.condition or "—",
        )
    console.print(table)
