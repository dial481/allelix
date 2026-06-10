# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Allelix command-line interface."""

from __future__ import annotations

import logging
import sys
import time
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple

import click
from rich.console import Console
from rich.table import Table

from allelix import __version__
from allelix.annotators import get_annotators
from allelix.databases import resolve_data_dir
from allelix.parsers import ParserNotFoundError, detect_parser, get_parser_by_name
from allelix.reports._pipeline import run_analysis
from allelix.reports.diff import compute_diff, load_previous_report
from allelix.reports.high_value import format_warnings, load_high_value_snps, scan_no_calls
from allelix.reports.html import render_html
from allelix.reports.json_report import render_json
from allelix.reports.methylation import METHYLATION_PANEL_GENES
from allelix.reports.terminal import render_terminal, render_terminal_diff

if TYPE_CHECKING:
    from allelix.annotators.base import Annotator
    from allelix.models import Variant
    from allelix.parsers.base import GenotypeParser

console = Console()

# Sort 1-22 numerically, then X, Y, MT, then anything else alphabetically.
_NAMED_CHROM_ORDER = {"X": 0, "Y": 1, "MT": 2}


def _chrom_sort_key(chrom: str) -> tuple[int, int, str]:
    """Sort key: autosomes (1-22), then X/Y/MT, then unknowns alphabetically."""
    if chrom.isdigit():
        return (0, int(chrom), "")
    if chrom in _NAMED_CHROM_ORDER:
        return (1, _NAMED_CHROM_ORDER[chrom], "")
    return (2, 0, chrom)


def _percent(part: int, total: int) -> str:
    if total == 0:
        return "0.00%"
    return f"{part / total * 100:.2f}%"


class _WarningCounter(logging.Handler):
    """Count warning records emitted by the parser pipeline."""

    def __init__(self) -> None:
        super().__init__(level=logging.WARNING)
        self.count = 0

    def emit(self, record: logging.LogRecord) -> None:
        self.count += 1


class _LoggerSnapshot(NamedTuple):
    """Captured state of a Python logger for restoration after CLI mutates it."""

    level: int
    propagate: bool


def _wire_parser_logging() -> tuple[_WarningCounter, logging.Handler, _LoggerSnapshot]:
    """Attach warning capture + stderr surfacing to the parsers logger."""
    parser_logger = logging.getLogger("allelix.parsers")
    counter = _WarningCounter()
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.WARNING)
    stderr_handler.setFormatter(logging.Formatter("warning: %(message)s"))
    snapshot = _LoggerSnapshot(level=parser_logger.level, propagate=parser_logger.propagate)
    parser_logger.addHandler(counter)
    parser_logger.addHandler(stderr_handler)
    parser_logger.setLevel(logging.WARNING)
    parser_logger.propagate = False
    return counter, stderr_handler, snapshot


def _unwire_parser_logging(
    counter: _WarningCounter,
    stderr_handler: logging.Handler,
    snapshot: _LoggerSnapshot,
) -> None:
    parser_logger = logging.getLogger("allelix.parsers")
    parser_logger.removeHandler(counter)
    parser_logger.removeHandler(stderr_handler)
    parser_logger.setLevel(snapshot.level)
    parser_logger.propagate = snapshot.propagate


def _resolve_parser(file_path: Path, fmt: str | None) -> GenotypeParser:
    try:
        return get_parser_by_name(fmt) if fmt else detect_parser(file_path)
    except ParserNotFoundError as exc:
        raise click.ClickException(str(exc)) from exc


def _ready_annotators(
    data_dir: Path | None,
    *,
    include_benign: bool = False,
    gwas_filter_traits: bool = True,
) -> tuple[Path, list[Annotator], list[Annotator]]:
    resolved = resolve_data_dir(data_dir)
    annotators = get_annotators(
        resolved, include_benign=include_benign, gwas_filter_traits=gwas_filter_traits
    )
    ready: list[Annotator] = []
    not_ready: list[Annotator] = []
    for a in annotators:
        if a.is_ready():
            ready.append(a)
        else:
            not_ready.append(a)
    if not ready:
        names = ", ".join(a.name for a in annotators)
        raise click.ClickException(
            f"No annotators are ready. Run `allelix db update` first. Registered: {names}"
        )
    return resolved, ready, not_ready


_STALENESS_SECONDS = 7 * 24 * 60 * 60  # 7 days


def _maybe_refresh_databases(data_dir: Path) -> None:
    """Check database mtimes; refresh any that are stale and have a changed remote signal.

    Only runs for annotators that download data (SNPedia excluded).
    If the network is unreachable, warns and continues with stale caches.
    """
    now = time.time()
    annotators = get_annotators(data_dir)
    for annotator in annotators:
        with annotator:
            if not annotator.requires_download or not annotator.is_ready():
                continue
            # Code-driven sources (commit-pinned HF caches) never change
            # at a fixed URL — skip the HEAD request. See ADR-0030.
            if not annotator.server_driven_freshness:
                continue
            db_files = list(data_dir.glob(f"{annotator.name}*sqlite*"))
            if not db_files:
                continue
            newest_mtime = max(f.stat().st_mtime for f in db_files)
            age = now - newest_mtime
            if age <= _STALENESS_SECONDS:
                continue

            remote = annotator.fetch_remote_signal()
            if remote is None:
                age_days = int(age / 86400)
                console.print(
                    f"[yellow]{annotator.display_name} database is {age_days} days old. "
                    "Run `allelix db update` when online.[/yellow]"
                )
                continue

            cached = annotator.cached_remote_signal()
            if cached == remote:
                continue

            console.print(f"[bold]Updating {annotator.display_name}…[/bold]")
            if _run_setup(annotator):
                console.print(
                    f"[green]✓ {annotator.display_name} updated[/green] "
                    f"(version {annotator.version() or '(unknown)'})"
                )


def _format_from_path(output: Path, override: str | None) -> str:
    if override:
        return override.lower()
    suffix = output.suffix.lower()
    if suffix == ".html":
        return "html"
    if suffix == ".json":
        return "json"
    raise click.ClickException(
        f"Cannot infer report format from {output.name!r}. "
        "Pass --report-format html|json explicitly."
    )


def _run_analysis_command(
    file_path: Path,
    fmt: str | None,
    data_dir: Path | None,
    output: Path | None,
    report_format: str | None,
    min_magnitude: float,
    category: str | None,
    genes: frozenset[str] | None,
    build: str | None = None,
    include_benign: bool = False,
    gwas_min_magnitude: float | None = None,
    snpedia_min_magnitude: float | None = None,
    exclude_sources: frozenset[str] | None = None,
    gwas_all: bool = False,
    diff_path: Path | None = None,
    no_update: bool = False,
    no_gnomad: bool = False,
    no_alphamissense: bool = False,
) -> None:
    resolved = resolve_data_dir(data_dir)
    if not no_update:
        _maybe_refresh_databases(resolved)
    parser = _resolve_parser(file_path, fmt)
    _, ready, not_ready = _ready_annotators(
        data_dir, include_benign=include_benign, gwas_filter_traits=not gwas_all
    )

    from allelix.config import load_config

    cfg = load_config(resolved)
    ready = [a for a in ready if cfg.is_enabled(a.name)]

    if exclude_sources:
        ready = [a for a in ready if a.name not in exclude_sources]

    gnomad_annotator = None
    if not no_gnomad:
        from allelix.annotators.gnomad import GnomadAnnotator

        for a in ready:
            if isinstance(a, GnomadAnnotator):
                gnomad_annotator = a
                break
    ready = [a for a in ready if a.name != "gnomad"]

    am_annotator = None
    if not no_alphamissense:
        from allelix.annotators.alphamissense import AlphaMissenseAnnotator

        for a in ready:
            if isinstance(a, AlphaMissenseAnnotator):
                am_annotator = a
                break
    ready = [a for a in ready if a.name != "alphamissense"]

    if not_ready:
        names = [a.name for a in not_ready]
        console.print(
            f"[yellow]Skipping unready annotators: {', '.join(names)}[/yellow] "
            "(run `allelix db update` to populate)"
        )

    all_active: list[Annotator] = list(ready)
    if gnomad_annotator is not None and gnomad_annotator.is_ready():
        all_active.append(gnomad_annotator)
    if am_annotator is not None and am_annotator.is_ready():
        all_active.append(am_annotator)
    versions = ", ".join(f"{a.display_name} ({a.version() or 'unknown'})" for a in all_active)
    console.print(f"[dim]Analyzing against: {versions}[/dim]")

    counter, stderr_handler, snapshot = _wire_parser_logging()
    try:
        result = run_analysis(
            file_path,
            parser,
            ready,
            skipped_count_provider=lambda: counter.count,
            build_override=build,
            gnomad=gnomad_annotator,
            alphamissense=am_annotator,
        )
    finally:
        _unwire_parser_logging(counter, stderr_handler, snapshot)

    _emit_build_diagnostics(result)

    high_value = load_high_value_snps()
    hv_rsids = set(high_value)
    hv_variants: list[Variant] = [v for v in parser.parse(file_path) if v.rsid in hv_rsids]
    hv_warnings = scan_no_calls(hv_variants, high_value)
    if hv_warnings:
        console.print(
            f"[bold red]Warning:[/bold red] {len(hv_warnings)} high-value SNP(s) returned no-call:"
        )
        for line in format_warnings(hv_warnings):
            console.print(f"  [red]⚠[/red] {line}")

    if counter.count:
        console.print(
            f"[yellow]Note:[/yellow] {counter.count:,} malformed line(s) skipped "
            "(see warnings on stderr)."
        )

    source_floors: dict[str, float] | None = None
    if gwas_min_magnitude is not None or snpedia_min_magnitude is not None:
        source_floors = {}
        if gwas_min_magnitude is not None:
            source_floors["gwas"] = gwas_min_magnitude
        if snpedia_min_magnitude is not None:
            source_floors["snpedia"] = snpedia_min_magnitude

    diff_result = None
    if diff_path is not None:
        try:
            prev = load_previous_report(diff_path)
        except ValueError as exc:
            raise click.ClickException(str(exc)) from exc
        filtered_for_diff = result.filter(
            min_magnitude=min_magnitude,
            category=category,
            genes=genes,
            source_min_magnitudes=source_floors,
        )
        from allelix.reports._pipeline import rollup_gwas_duplicates

        filtered_for_diff = rollup_gwas_duplicates(filtered_for_diff)
        diff_result = compute_diff(
            filtered_for_diff,
            prev["annotations"],
            prev.get("generated_at", ""),
        )

    if output is None:
        if diff_result is not None:
            rendered = render_terminal_diff(diff_result, console)
        else:
            rendered = render_terminal(
                result,
                console=console,
                min_magnitude=min_magnitude,
                category=category,
                genes=genes,
                source_min_magnitudes=source_floors,
            )
    else:
        chosen = _format_from_path(output, report_format)
        hv_warning_lines = format_warnings(hv_warnings) if hv_warnings else None
        if chosen == "json":
            hv_dicts = (
                [{"rsid": w.snp.rsid, "gene": w.snp.gene, "note": w.snp.note} for w in hv_warnings]
                if hv_warnings
                else None
            )
            rendered = render_json(
                result,
                output_path=output,
                min_magnitude=min_magnitude,
                category=category,
                genes=genes,
                source_min_magnitudes=source_floors,
                diff=diff_result,
                high_value_no_calls=hv_dicts,
            )
        else:
            rendered = render_html(
                result,
                output_path=output,
                min_magnitude=min_magnitude,
                category=category,
                genes=genes,
                source_min_magnitudes=source_floors,
                diff=diff_result,
                high_value_no_calls=hv_warning_lines,
            )
        console.print(f"[green]Wrote {rendered:,} annotation(s) to {output}[/green]")

    console.print(
        f"[dim]{len(result.annotations):,} total annotation(s) from {len(ready)} "
        f"database(s) across {result.total_variants:,} variant(s).[/dim]"
    )


@click.group()
@click.version_option(version=__version__, prog_name="allelix")
def main() -> None:
    """Allelix: open-source genotype analysis toolkit."""


@main.command()
@click.argument(
    "file_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--format",
    "fmt",
    default=None,
    help="Force a specific parser (e.g., myhappygenes). Default: auto-detect.",
)
def stats(file_path: Path, fmt: str | None) -> None:
    """Show summary statistics for a genotype file."""
    parser = _resolve_parser(file_path, fmt)
    counter, stderr_handler, snapshot = _wire_parser_logging()

    high_value = load_high_value_snps()
    hv_rsids = set(high_value)
    hv_variants: list[Variant] = []

    total = 0
    no_calls = 0
    het = 0
    hom = 0
    chrom_counts: dict[str, int] = {}
    try:
        metadata = parser.get_metadata(file_path)
        for variant in parser.parse(file_path):
            total += 1
            if variant.rsid in hv_rsids:
                hv_variants.append(variant)
            if variant.is_no_call:
                no_calls += 1
            elif variant.is_heterozygous:
                het += 1
            else:
                hom += 1
            chrom_counts[variant.chromosome] = chrom_counts.get(variant.chromosome, 0) + 1
    finally:
        _unwire_parser_logging(counter, stderr_handler, snapshot)

    summary = Table(title=f"Genotype File Stats: {file_path.name}")
    summary.add_column("Metric", style="cyan", no_wrap=True)
    summary.add_column("Value", style="green")
    summary.add_row("Format", parser.display_name)
    summary.add_row("Sample ID", metadata["sample_id"] or "(unknown)")
    summary.add_row("Build", metadata["build"])
    summary.add_row("Total SNPs", f"{total:,}")
    summary.add_row("No-calls", f"{no_calls:,} ({_percent(no_calls, total)})")
    summary.add_row("Heterozygous", f"{het:,} ({_percent(het, total)})")
    summary.add_row("Homozygous", f"{hom:,} ({_percent(hom, total)})")
    if counter.count:
        summary.add_row(
            "Skipped (malformed)",
            f"[yellow]{counter.count:,}[/yellow] (see warnings on stderr)",
        )

    hv_warnings = scan_no_calls(hv_variants, high_value)
    if hv_warnings:
        summary.add_row(
            "High-value no-calls",
            f"[red]{len(hv_warnings)}[/red]",
        )
    console.print(summary)

    if hv_warnings:
        for line in format_warnings(hv_warnings):
            console.print(f"  [red]⚠[/red] {line}")

    chrom_table = Table(title="Variants per Chromosome")
    chrom_table.add_column("Chromosome", style="cyan", no_wrap=True)
    chrom_table.add_column("Count", style="green", justify="right")
    for chrom in sorted(chrom_counts, key=_chrom_sort_key):
        chrom_table.add_row(chrom, f"{chrom_counts[chrom]:,}")
    console.print(chrom_table)


_FILE_ARG = click.argument(
    "file_path", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
_FORMAT_OPT = click.option(
    "--format", "fmt", default=None, help="Force a specific parser. Default: auto-detect."
)
_DATA_DIR_OPT = click.option(
    "--data-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Override database cache location.",
)
_MIN_MAG_OPT = click.option(
    "--min-magnitude",
    type=float,
    default=5.0,
    show_default=True,
    help="Filter annotations below this magnitude. Use 0 for the full unfiltered set.",
)
_OUTPUT_OPT = click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Write a report file (.html or .json). Omit for terminal output.",
)
_REPORT_FORMAT_OPT = click.option(
    "--report-format",
    type=click.Choice(["html", "json"], case_sensitive=False),
    default=None,
    help="Override report format detection (otherwise inferred from --output extension).",
)
_INCLUDE_BENIGN_OPT = click.option(
    "--include-benign",
    is_flag=True,
    default=False,
    help="Include ClinVar Benign/Likely_benign annotations (suppressed by default).",
)
_GWAS_MIN_MAG_OPT = click.option(
    "--gwas-min-magnitude",
    type=float,
    default=9.0,
    show_default=True,
    help="Magnitude floor for GWAS Catalog annotations (overrides --min-magnitude for GWAS).",
)
_SNPEDIA_MIN_MAG_OPT = click.option(
    "--snpedia-min-magnitude",
    type=float,
    default=2.0,
    show_default=True,
    help="Magnitude floor for SNPedia annotations (overrides --min-magnitude for SNPedia).",
)
_INCLUDE_GWAS_OPT = click.option(
    "--include-gwas",
    is_flag=True,
    default=False,
    help="Include GWAS Catalog annotations (excluded by default in focused reports).",
)
_EXCLUDE_SNPEDIA_OPT = click.option(
    "--exclude-snpedia",
    is_flag=True,
    default=False,
    help="Exclude SNPedia annotations. Required for commercial use (CC BY-NC-SA 3.0).",
)
_GWAS_ALL_OPT = click.option(
    "--gwas-all",
    is_flag=True,
    default=False,
    help="Include all GWAS trait categories (disables default noise filtering).",
)
_DIFF_OPT = click.option(
    "--diff",
    "diff_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Dev/QA tool: compare current output against a previous JSON report "
        "to detect regressions from code changes, database refreshes, or "
        "filter adjustments. Shows new, changed, and removed annotations. "
        "Not a monitoring tool — use for version-to-version validation."
    ),
)
_NO_UPDATE_OPT = click.option(
    "--no-update",
    is_flag=True,
    default=False,
    help="Skip the pre-analysis database freshness check.",
)
_NO_GNOMAD_OPT = click.option(
    "--no-gnomad",
    is_flag=True,
    default=False,
    help="Skip gnomAD population frequency enrichment.",
)
_NO_ALPHAMISSENSE_OPT = click.option(
    "--no-alphamissense",
    is_flag=True,
    default=False,
    help="Skip AlphaMissense variant pathogenicity enrichment.",
)
_BUILD_OPT = click.option(
    "--build",
    type=click.Choice(["grch37", "grch38", "auto"], case_sensitive=False),
    default="auto",
    help=(
        "Genome build of the input file. 'auto' detects from position data "
        "(ADR-0021) and ignores the file header. 'grch37' / 'grch38' force a "
        "specific build, skipping detection."
    ),
)


def _resolve_clinvar_builds(value: str) -> tuple[str, ...]:
    """Map a `db update --build` value to a tuple of build identifiers."""
    v = (value or "both").strip().lower()
    if v == "both":
        return ("GRCh37", "GRCh38")
    if v == "grch37":
        return ("GRCh37",)
    if v == "grch38":
        return ("GRCh38",)
    raise click.ClickException(f"Unknown --build value {value!r}")


def _normalize_cli_build(value: str | None) -> str | None:
    """Map a --build CLI value to a canonical build identifier or None for auto."""
    if value is None:
        return None
    v = value.strip().lower()
    if v in ("", "auto"):
        return None
    if v == "grch37":
        return "GRCh37"
    if v == "grch38":
        return "GRCh38"
    raise click.ClickException(f"Unknown --build value {value!r}")


def _emit_build_diagnostics(result: object) -> None:
    """Print a one-line build banner and a warning on header/data mismatch."""
    diag = getattr(result, "build_diagnostics", None)
    if diag is None:
        return
    matched = f"{diag.matched_count}/{diag.inspected_count}" if diag.inspected_count else "0/0"
    if diag.override:
        source = "override"
    elif diag.detected_build:
        source = "detected"
    elif diag.header_build:
        source = "header (no position confirmation)"
    else:
        source = "fallback (no known SNPs matched)"
    console.print(
        f"[dim]Build: {diag.effective_build} ({source}; "
        f"{matched} known-SNP positions matched)[/dim]"
    )
    if diag.mismatch:
        console.print(
            f"[yellow]Build mismatch: file header claims {diag.header_build} but "
            f"position data is {diag.detected_build}. Using {diag.detected_build}. "
            f"This is a real-world data-quality issue — your provider may have "
            f"mislabeled the build (see ADR-0021).[/yellow]"
        )
    if diag.effective_build == "GRCh36":
        console.print(
            "[yellow]Warning: GRCh36 (hg18) detected. rsID-based annotations "
            "(PharmGKB, GWAS Catalog, SNPedia, gnomAD) are complete. ClinVar "
            "position-matching is skipped (no GRCh36 cache — see ADR-0025). "
            "For full ClinVar coverage, liftOver to GRCh38 first: "
            "docs/grch36-liftover.md[/yellow]"
        )


@main.command()
@_FILE_ARG
@_FORMAT_OPT
@_DATA_DIR_OPT
@_MIN_MAG_OPT
@click.option(
    "--category",
    type=str,
    default=None,
    help="Filter to a single bucket (clinical, pharma).",
)
@_OUTPUT_OPT
@_REPORT_FORMAT_OPT
@_BUILD_OPT
@_INCLUDE_BENIGN_OPT
@_GWAS_MIN_MAG_OPT
@_SNPEDIA_MIN_MAG_OPT
@_GWAS_ALL_OPT
@_EXCLUDE_SNPEDIA_OPT
@_DIFF_OPT
@_NO_UPDATE_OPT
@_NO_GNOMAD_OPT
@_NO_ALPHAMISSENSE_OPT
def analyze(
    file_path: Path,
    fmt: str | None,
    data_dir: Path | None,
    min_magnitude: float,
    category: str | None,
    output: Path | None,
    report_format: str | None,
    build: str,
    include_benign: bool,
    gwas_min_magnitude: float,
    snpedia_min_magnitude: float,
    gwas_all: bool,
    exclude_snpedia: bool,
    diff_path: Path | None,
    no_update: bool,
    no_gnomad: bool,
    no_alphamissense: bool,
) -> None:
    """Annotate a genotype file against all ready reference databases."""
    _run_analysis_command(
        file_path=file_path,
        fmt=fmt,
        data_dir=data_dir,
        output=output,
        report_format=report_format,
        min_magnitude=min_magnitude,
        category=category,
        genes=None,
        build=_normalize_cli_build(build),
        include_benign=include_benign,
        gwas_min_magnitude=gwas_min_magnitude,
        snpedia_min_magnitude=snpedia_min_magnitude,
        exclude_sources=frozenset({"snpedia"}) if exclude_snpedia else None,
        gwas_all=gwas_all,
        diff_path=diff_path,
        no_update=no_update,
        no_gnomad=no_gnomad,
        no_alphamissense=no_alphamissense,
    )


@main.command()
@_FILE_ARG
@_FORMAT_OPT
@click.option(
    "--snps",
    required=True,
    help="Comma-separated rsIDs to extract (e.g., rs1801133,rs4680).",
)
def extract(file_path: Path, fmt: str | None, snps: str) -> None:
    """Print diploid genotypes for specific rsIDs — spot-check carrier status.

    Useful for verifying ClinVar / PharmGKB hits against the actual file
    before trusting them. The "Genotype" column shows the diploid call as
    the array (or VCF) reported it; "Het?" and "No-call?" answer the
    questions the carrier rule (ADR-0007) actually checks.
    """
    parser = _resolve_parser(file_path, fmt)
    wanted = {s.strip() for s in snps.split(",") if s.strip()}
    if not wanted:
        raise click.ClickException("--snps cannot be empty.")

    counter, stderr_handler, snapshot = _wire_parser_logging()
    found: dict[str, object] = {}
    try:
        for variant in parser.parse(file_path):
            if variant.rsid in wanted:
                found[variant.rsid] = variant
                if len(found) == len(wanted):
                    break  # streaming early-exit once we have everything
    finally:
        _unwire_parser_logging(counter, stderr_handler, snapshot)

    table = Table(title=f"Genotypes from {file_path.name}")
    table.add_column("rsID", style="cyan", no_wrap=True)
    table.add_column("Chr", no_wrap=True)
    table.add_column("Position", justify="right")
    table.add_column("Genotype", style="yellow", no_wrap=True)
    table.add_column("Het?", justify="center")
    table.add_column("No-call?", justify="center")
    for rsid in sorted(wanted):
        variant = found.get(rsid)
        if variant is None:
            table.add_row(rsid, "—", "—", "[red]not in file[/red]", "—", "—")
            continue
        table.add_row(
            variant.rsid,
            variant.chromosome,
            f"{variant.position:,}",
            variant.genotype,
            "yes" if variant.is_heterozygous else "no",
            "[red]yes[/red]" if variant.is_no_call else "no",
        )
    console.print(table)


@main.command()
@click.argument("file1", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("file2", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--format1", "fmt1", default=None, help="Force parser for file 1.")
@click.option("--format2", "fmt2", default=None, help="Force parser for file 2.")
def compare(file1: Path, file2: Path, fmt1: str | None, fmt2: str | None) -> None:
    """Compare two genotype files — coverage overlap and concordance.

    Reports shared rsIDs, file-specific rsIDs, genotype agreement,
    strand-flip matches (complementary alleles on opposite strands),
    discordant calls, and strand-ambiguous positions.
    """
    from allelix.compare import compare_variants
    from allelix.utils.build_detect import detect_build

    parser1 = _resolve_parser(file1, fmt1)
    parser2 = _resolve_parser(file2, fmt2)
    variants1 = list(parser1.parse(file1))
    variants2 = list(parser2.parse(file2))

    det1 = detect_build(variants1)
    det2 = detect_build(variants2)
    build1 = det1.build or parser1.get_metadata(file1).get("build", "unknown")
    build2 = det2.build or parser2.get_metadata(file2).get("build", "unknown")

    result = compare_variants(variants1, variants2, build1=build1, build2=build2)

    if result.build1 != result.build2:
        console.print(
            f"[yellow]Warning: builds differ ({result.build1} vs {result.build2}). "
            "Position-based comparisons may be unreliable.[/yellow]"
        )

    table = Table(title="Coverage Summary")
    table.add_column("Metric", style="bold")
    table.add_column("Value", justify="right")
    table.add_row("File 1", f"{file1.name} ({result.file1_total:,} variants)")
    table.add_row("File 2", f"{file2.name} ({result.file2_total:,} variants)")
    table.add_row("Build (file 1)", result.build1)
    table.add_row("Build (file 2)", result.build2)
    table.add_row("Shared rsIDs", f"{result.shared:,}")
    table.add_row("File 1 only", f"{result.file1_only:,}")
    table.add_row("File 2 only", f"{result.file2_only:,}")
    console.print(table)

    conc_table = Table(title="Genotype Concordance")
    conc_table.add_column("Category", style="bold")
    conc_table.add_column("Count", justify="right")
    conc_table.add_column("%", justify="right")
    for label, count in [
        ("Concordant", result.concordant),
        ("Strand-flip match", result.strand_flip_match),
        ("Discordant", result.discordant),
        ("Strand-ambiguous", result.strand_ambiguous),
        ("No-call (either file)", result.no_call),
    ]:
        pct = _percent(count, result.shared) if result.shared else "—"
        conc_table.add_row(label, f"{count:,}", pct)
    console.print(conc_table)

    if result.chromosome_counts:
        chrom_table = Table(title="Per-Chromosome Breakdown")
        chrom_table.add_column("Chr", style="cyan", no_wrap=True)
        chrom_table.add_column("Concordant", justify="right")
        chrom_table.add_column("Flip", justify="right")
        chrom_table.add_column("Discordant", justify="right")
        chrom_table.add_column("Ambiguous", justify="right")
        chrom_table.add_column("No-call", justify="right")
        for chrom in sorted(result.chromosome_counts, key=_chrom_sort_key):
            c = result.chromosome_counts[chrom]
            chrom_table.add_row(
                chrom,
                str(c.get("concordant", 0)),
                str(c.get("strand_flip_match", 0)),
                str(c.get("discordant", 0)),
                str(c.get("strand_ambiguous", 0)),
                str(c.get("no_call", 0)),
            )
        console.print(chrom_table)


@main.command()
@_FILE_ARG
@_FORMAT_OPT
@_DATA_DIR_OPT
@_MIN_MAG_OPT
@_OUTPUT_OPT
@_REPORT_FORMAT_OPT
@_BUILD_OPT
@_INCLUDE_BENIGN_OPT
@_GWAS_MIN_MAG_OPT
@_SNPEDIA_MIN_MAG_OPT
@_INCLUDE_GWAS_OPT
@_GWAS_ALL_OPT
@_EXCLUDE_SNPEDIA_OPT
@_DIFF_OPT
@_NO_UPDATE_OPT
@_NO_GNOMAD_OPT
@_NO_ALPHAMISSENSE_OPT
def methylation(
    file_path: Path,
    fmt: str | None,
    data_dir: Path | None,
    min_magnitude: float,
    output: Path | None,
    report_format: str | None,
    build: str,
    include_benign: bool,
    gwas_min_magnitude: float,
    snpedia_min_magnitude: float,
    include_gwas: bool,
    gwas_all: bool,
    exclude_snpedia: bool,
    diff_path: Path | None,
    no_update: bool,
    no_gnomad: bool,
    no_alphamissense: bool,
) -> None:
    """Methylation-pathway-focused report (MTHFR, MTR, MTRR, COMT, CBS, …)."""
    excluded: set[str] = set()
    if not include_gwas:
        excluded.add("gwas")
    if exclude_snpedia:
        excluded.add("snpedia")
    _run_analysis_command(
        file_path=file_path,
        fmt=fmt,
        data_dir=data_dir,
        output=output,
        report_format=report_format,
        min_magnitude=min_magnitude,
        category=None,
        genes=METHYLATION_PANEL_GENES,
        build=_normalize_cli_build(build),
        include_benign=include_benign,
        gwas_min_magnitude=gwas_min_magnitude,
        snpedia_min_magnitude=snpedia_min_magnitude,
        exclude_sources=frozenset(excluded) if excluded else None,
        gwas_all=gwas_all,
        diff_path=diff_path,
        no_update=no_update,
        no_gnomad=no_gnomad,
        no_alphamissense=no_alphamissense,
    )


@main.command()
@_FILE_ARG
@_FORMAT_OPT
@_DATA_DIR_OPT
@_MIN_MAG_OPT
@_OUTPUT_OPT
@_REPORT_FORMAT_OPT
@_BUILD_OPT
@_INCLUDE_BENIGN_OPT
@_GWAS_MIN_MAG_OPT
@_SNPEDIA_MIN_MAG_OPT
@_INCLUDE_GWAS_OPT
@_GWAS_ALL_OPT
@_EXCLUDE_SNPEDIA_OPT
@_DIFF_OPT
@_NO_UPDATE_OPT
@_NO_GNOMAD_OPT
@_NO_ALPHAMISSENSE_OPT
def pharmacogenomics(
    file_path: Path,
    fmt: str | None,
    data_dir: Path | None,
    min_magnitude: float,
    output: Path | None,
    report_format: str | None,
    build: str,
    include_benign: bool,
    gwas_min_magnitude: float,
    snpedia_min_magnitude: float,
    include_gwas: bool,
    gwas_all: bool,
    exclude_snpedia: bool,
    diff_path: Path | None,
    no_update: bool,
    no_gnomad: bool,
    no_alphamissense: bool,
) -> None:
    """Pharmacogenomics-focused report (annotations from PharmGKB-style sources)."""
    excluded: set[str] = set()
    if not include_gwas:
        excluded.add("gwas")
    if exclude_snpedia:
        excluded.add("snpedia")
    _run_analysis_command(
        file_path=file_path,
        fmt=fmt,
        data_dir=data_dir,
        output=output,
        report_format=report_format,
        min_magnitude=min_magnitude,
        category="pharma",
        genes=None,
        build=_normalize_cli_build(build),
        include_benign=include_benign,
        gwas_min_magnitude=gwas_min_magnitude,
        snpedia_min_magnitude=snpedia_min_magnitude,
        exclude_sources=frozenset(excluded) if excluded else None,
        gwas_all=gwas_all,
        diff_path=diff_path,
        no_update=no_update,
        no_gnomad=no_gnomad,
        no_alphamissense=no_alphamissense,
    )


@main.group()
def db() -> None:
    """Manage local reference database cache."""


def _stamp_remote_signal(annotator: Annotator, signal: str) -> None:
    """Write a remote signal to an existing cache without re-downloading."""
    import contextlib
    import sqlite3

    from allelix.databases.manager import stamp_remote_signal

    db_path = getattr(annotator, "_db_path", None)
    if db_path is None:
        return
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        stamp_remote_signal(conn, annotator.name, signal)
        conn.commit()


def _run_setup(annotator: Annotator) -> bool:
    """Invoke annotator.setup(). Returns True on success, False on failure."""
    try:
        annotator.setup()
    except Exception as exc:
        if hasattr(exc, "close"):
            exc.close()
        console.print(f"  [red]{annotator.name}: {exc}[/red]")
        return False
    sig = getattr(annotator, "cached_remote_signal", lambda: None)()
    if sig and "cpic:unavailable" in sig:
        console.print(
            f"  [yellow]{annotator.name}: updated (CPIC unavailable — "
            "non-finding filter degraded, retry later)[/yellow]"
        )
    return True


@db.command("update")
@_DATA_DIR_OPT
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Re-download even if the local cache appears current.",
)
@click.option(
    "--no-gnomad",
    is_flag=True,
    default=False,
    help="Skip gnomAD population frequency database.",
)
@click.option(
    "--no-alphamissense",
    is_flag=True,
    default=False,
    help="Skip AlphaMissense pathogenicity database.",
)
@click.option(
    "--build",
    type=click.Choice(["grch37", "grch38", "both"], case_sensitive=False),
    default="both",
    help=(
        "Which ClinVar genome build(s) to download. 'both' (default) keeps "
        "GRCh37 and GRCh38 caches in sync so `analyze` can dispatch by "
        "detected build (ADR-0021). 'grch37' / 'grch38' restrict to one to "
        "save bandwidth."
    ),
)
def db_update(
    data_dir: Path | None, force: bool, no_gnomad: bool, no_alphamissense: bool, build: str
) -> None:
    """Download or refresh reference databases.

    For each annotator:
      - no cache → download
      - --force → download
      - cache + remote signal matches cache → skip
      - cache + remote signal differs (or legacy v0.4.1 cache with no
        stored signal) → download
      - cache + remote signal can't be fetched → skip with notice (use
        --force to override)

    `--build` selects which ClinVar build(s) to manage. Default 'both'
    downloads GRCh37 and GRCh38 caches.
    """
    resolved = resolve_data_dir(data_dir)
    console.print(f"Data directory: [cyan]{resolved}[/cyan]")
    clinvar_builds = _resolve_clinvar_builds(build)
    for annotator in get_annotators(resolved, clinvar_builds=clinvar_builds):
        with annotator:
            if no_gnomad and annotator.name == "gnomad":
                console.print(f"  [dim]{annotator.name}: skipped (--no-gnomad)[/dim]")
                continue
            if no_alphamissense and annotator.name == "alphamissense":
                console.print(f"  [dim]{annotator.name}: skipped (--no-alphamissense)[/dim]")
                continue

            if not annotator.requires_download:
                if annotator.is_ready():
                    console.print(
                        f"  [dim]{annotator.name}: ready "
                        f"({annotator.version() or 'unknown'})[/dim]"
                    )
                continue

            if not annotator.is_ready():
                console.print(f"  [bold]{annotator.name}[/bold]: downloading…")
                if _run_setup(annotator):
                    console.print(
                        f"  [green]✓ {annotator.name} ready[/green] "
                        f"(version {annotator.version() or '(unknown)'})"
                    )
                continue

            if force:
                console.print(f"  [bold]{annotator.name}[/bold]: --force; refreshing…")
                if _run_setup(annotator):
                    console.print(
                        f"  [green]✓ {annotator.name} refreshed[/green] "
                        f"(version {annotator.version() or '(unknown)'})"
                    )
                continue

            # Code-driven sources (commit-pinned HF caches) are updated
            # only via code changes — no runtime freshness probe needed.
            if not annotator.server_driven_freshness:
                console.print(
                    f"  [dim]{annotator.name}: already current "
                    f"(version {annotator.version() or '(unknown)'})[/dim]"
                )
                continue

            remote = annotator.fetch_remote_signal()
            if remote is None:
                console.print(
                    f"  [yellow]{annotator.name}: cache present, but remote "
                    "freshness can't be verified (network error or no signal). "
                    "Pass --force to refresh anyway.[/yellow]"
                )
                continue

            cached = annotator.cached_remote_signal()
            if cached == remote:
                console.print(
                    f"  [dim]{annotator.name}: already current "
                    f"(version {annotator.version() or '(unknown)'})[/dim]"
                )
                continue

            if cached is None:
                _stamp_remote_signal(annotator, remote)
                console.print(
                    f"  [dim]{annotator.name}: stamped remote signal "
                    f"(version {annotator.version() or '(unknown)'})[/dim]"
                )
                continue

            console.print(f"  [bold]{annotator.name}[/bold]: remote signal changed; refreshing…")
            if _run_setup(annotator):
                console.print(
                    f"  [green]✓ {annotator.name} refreshed[/green] "
                    f"(version {annotator.version() or '(unknown)'})"
                )


@db.command("status")
@_DATA_DIR_OPT
def db_status(data_dir: Path | None) -> None:
    """Show installed reference database versions and freshness."""
    resolved = resolve_data_dir(data_dir)
    table = Table(title=f"Reference Databases ({resolved})")
    table.add_column("Annotator", style="cyan", no_wrap=True)
    table.add_column("Ready", justify="center")
    table.add_column("Version")
    table.add_column("Records", justify="right")
    for annotator in get_annotators(resolved):
        with annotator:
            ready = annotator.is_ready()
            ready_marker = "[green]yes[/green]" if ready else "[red]no[/red]"
            version = annotator.version() or "—"
            sig = getattr(annotator, "cached_remote_signal", lambda: None)()
            if sig and "cpic:unavailable" in sig:
                version += " (no CPIC)"
            records = "—"
            count_fn = getattr(annotator, "record_count", None)
            if callable(count_fn):
                count = count_fn()
                if count is not None:
                    records = f"{count:,}"
            table.add_row(annotator.display_name, ready_marker, version, records)
    console.print(table)


@main.group()
def config() -> None:
    """Manage persistent configuration (source toggles, license mode)."""


@config.command("show")
@_DATA_DIR_OPT
def config_show(data_dir: Path | None) -> None:
    """Display current configuration."""
    from allelix.config import NON_COMMERCIAL_SOURCES, load_config

    resolved = resolve_data_dir(data_dir)
    cfg = load_config(resolved)

    table = Table(title=f"Configuration ({resolved / 'config.toml'})")
    table.add_column("Source", style="cyan", no_wrap=True)
    table.add_column("Enabled", justify="center")
    table.add_column("Note", style="dim")
    for name, enabled in sorted(cfg.sources.items()):
        if cfg.commercial and name in NON_COMMERCIAL_SOURCES:
            marker = "[red]no[/red]"
            note = "disabled by commercial mode"
        elif enabled:
            marker = "[green]yes[/green]"
            note = ""
        else:
            marker = "[red]no[/red]"
            note = ""
        table.add_row(name, marker, note)
    console.print(table)
    mode = "[yellow]commercial[/yellow]" if cfg.commercial else "[green]personal[/green]"
    console.print(f"License mode: {mode}")


@config.command("set")
@_DATA_DIR_OPT
@click.argument("key")
@click.argument("value")
def config_set(data_dir: Path | None, key: str, value: str) -> None:
    r"""Set a configuration value.

    \b
    Keys:
      sources.<name>     Enable/disable a source (true/false)
      license.commercial Set commercial mode (true/false)

    \b
    Examples:
      allelix config set sources.snpedia false
      allelix config set license.commercial true
    """
    from allelix.config import load_config, save_config

    resolved = resolve_data_dir(data_dir)
    cfg = load_config(resolved)

    val_lower = value.strip().lower()
    if val_lower not in ("true", "false"):
        raise click.ClickException(f"Value must be 'true' or 'false', got {value!r}")
    bool_val = val_lower == "true"

    if key.startswith("sources."):
        source_name = key[len("sources.") :]
        cfg.sources[source_name] = bool_val
    elif key == "license.commercial":
        cfg.commercial = bool_val
    else:
        raise click.ClickException(
            f"Unknown key {key!r}. Use 'sources.<name>' or 'license.commercial'."
        )

    save_config(resolved, cfg)
    console.print(f"[green]Set {key} = {val_lower}[/green]")


if __name__ == "__main__":
    main()
