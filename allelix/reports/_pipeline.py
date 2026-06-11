# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Shared analysis pipeline used by `analyze`, `methylation`, and `pharmacogenomics`.

The CLI builds an `AnalysisResult` once and hands it to a renderer
(terminal, JSON, HTML). Renderers never query the database or re-iterate
the parser — they receive a fully-populated value object.

ADR-0021: this pipeline owns build detection. Parsers report the
header-claimed build; the pipeline replaces each variant's `build`
with the build detected from position data (or the user's `--build`
override) before annotators see the variant.
"""

from __future__ import annotations

import contextlib
import logging
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

from allelix.utils.build_detect import (
    BUILD_GRCH36,
    BUILD_GRCH37,
    BUILD_GRCH38,
    KNOWN_SNP_POSITIONS,
    detect_build,
    normalize_build_label,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable
    from pathlib import Path

    from allelix.annotators.alphamissense import AlphaMissenseAnnotator
    from allelix.annotators.base import Annotator
    from allelix.annotators.cadd import CaddAnnotator
    from allelix.annotators.gnomad import GnomadAnnotator
    from allelix.models import Annotation, Variant
    from allelix.parsers.base import GenotypeParser


# How many input variants to buffer while waiting for detection to
# converge. Detection completes once every entry in KNOWN_SNP_POSITIONS
# has been seen; typical files cover the table within the first ~5000
# probes. Cap so a file with no known SNPs doesn't buffer the whole
# input.
_DETECTION_BUFFER_LIMIT = 100_000


@dataclass
class BuildDiagnostics:
    """What the pipeline learned about the file's genome build.

    `header_build` is the build claimed by the file header (normalized
    to GRCh37/GRCh38 via `normalize_build_label`; may be None if the
    header doesn't say or uses an unrecognized label).

    `detected_build` is what position data says (None if no known SNPs
    appeared in the input).

    `effective_build` is what was actually used for annotation — either
    a CLI `--build` override, the detected build, or a fallback. Always
    set when the pipeline ran on any data.

    `mismatch` is True when header_build and detected_build disagree
    AND no override was supplied. The CLI surfaces this as a warning.
    """

    header_build: str | None
    detected_build: str | None
    effective_build: str
    override: bool
    matched_count: int
    inspected_count: int

    @property
    def mismatch(self) -> bool:
        return (
            not self.override
            and self.header_build is not None
            and self.detected_build is not None
            and self.header_build != self.detected_build
        )


@dataclass
class AnalysisResult:
    """Everything a renderer needs to produce a report."""

    file_path: Path
    parser_name: str
    parser_display_name: str
    sample_id: str
    build: str
    total_variants: int
    skipped_count: int
    annotators_used: list[tuple[str, str | None]]
    annotations: list[Annotation] = field(default_factory=list)
    build_diagnostics: BuildDiagnostics | None = None

    def filter(
        self,
        *,
        min_magnitude: float = 0.0,
        category: str | None = None,
        genes: Iterable[str] | None = None,
        source_min_magnitudes: dict[str, float] | None = None,
    ) -> list[Annotation]:
        """Apply the standard filters and return a sorted list of annotations.

        Filters are independent and combine with AND. Sort is by magnitude
        descending, then rsid ascending (stable, deterministic).

        `source_min_magnitudes` overrides the floor for specific sources
        (e.g. ``{"gwas": 9.0, "snpedia": 2.0}``). When a source has an
        entry, that value IS the floor for that source — it can raise OR
        lower the global ``min_magnitude``. Sources without an entry use
        the global floor.
        """
        gene_set = {g.upper() for g in genes} if genes else None
        out: list[Annotation] = []
        for a in self.annotations:
            if (
                source_min_magnitudes
                and a.source in source_min_magnitudes
                and not a.is_must_include
            ):
                floor = source_min_magnitudes[a.source]
            else:
                floor = min_magnitude
            if a.magnitude < floor:
                continue
            if category is not None and a.category != category:
                continue
            if gene_set is not None and (a.gene or "").upper() not in gene_set:
                continue
            out.append(a)
        out.sort(key=lambda a: (-a.magnitude, a.rsid))
        return out


def _gwas_base_trait(description: str) -> str | None:
    """Extract trait text from a GWAS description, stripping MTAG suffix and PheCode label."""
    marker = "GWAS Catalog: "
    if marker not in description:
        return None
    s = description.split(marker, 1)[1]
    s = s.split(" (p=", 1)[0]
    if s.endswith(" (MTAG)"):
        s = s[: -len(" (MTAG)")]
    s = s.split(" (PheCode ", 1)[0]
    return s.strip().lower()


def _gwas_phecode_parent(description: str) -> str | None:
    """Extract PheCode parent (numeric prefix before the dot), or None."""
    idx = description.find("(PheCode ")
    if idx == -1:
        return None
    rest = description[idx + len("(PheCode ") :]
    end = rest.find(")")
    if end == -1:
        return None
    code = rest[:end].strip()
    parent = code.split(".", 1)[0]
    return parent if parent.isdigit() else None


def _gwas_p_value(description: str) -> float:
    """Extract p-value from a GWAS description. Returns inf if unparseable."""
    idx = description.find("(p=")
    if idx == -1:
        return float("inf")
    rest = description[idx + len("(p=") :]
    end = rest.find(",")
    if end == -1:
        end = rest.find(")")
    if end == -1:
        return float("inf")
    try:
        return float(rest[:end].strip())
    except ValueError:
        return float("inf")


def rollup_gwas_duplicates(annotations: list[Annotation]) -> list[Annotation]:
    """Collapse GWAS MTAG twins and PheCode parent/child hierarchies.

    Operates on the filtered annotation list (the output of
    AnalysisResult.filter). Non-GWAS rows pass through untouched.
    Must-include rows are never dropped.

    See ADR-0024 'MTAG and PheCode rollup' for rules.
    """
    survivors: list[Annotation] = []
    gwas_rows: list[Annotation] = []
    for a in annotations:
        (gwas_rows if a.source == "gwas" else survivors).append(a)

    if not gwas_rows:
        return annotations

    plain_keys = {
        (a.rsid, _gwas_base_trait(a.description))
        for a in gwas_rows
        if "(MTAG)" not in a.description
    }
    after_mtag = [
        a
        for a in gwas_rows
        if a.is_must_include
        or "(MTAG)" not in a.description
        or (a.rsid, _gwas_base_trait(a.description)) not in plain_keys
    ]

    by_parent: dict[tuple[str, str], list[Annotation]] = {}
    no_phecode: list[Annotation] = []
    for a in after_mtag:
        parent = _gwas_phecode_parent(a.description)
        if parent is None or a.is_must_include:
            no_phecode.append(a)
        else:
            by_parent.setdefault((a.rsid, parent), []).append(a)
    for group in by_parent.values():
        winner = min(group, key=lambda x: _gwas_p_value(x.description))
        no_phecode.append(winner)

    survivors.extend(no_phecode)
    survivors.sort(key=lambda a: (-a.magnitude, a.rsid))
    return survivors


def _enrich_cadd(
    annotations: list[Annotation],
    gnomad: GnomadAnnotator,
    cadd: CaddAnnotator,
) -> None:
    """Stamp annotations with CADD PHRED scores via coordinate resolution.

    Resolves rsIDs to genomic coordinates through gnomAD, normalizes
    alleles to reference-forward orientation, and looks up CADD scores.
    """
    from allelix.utils.allele import resolve_strand

    rsids = {a.rsid for a in annotations}
    coord_map = gnomad.bulk_resolve_coordinates(rsids)
    if not coord_map:
        return

    cadd_keys: set[tuple[str, int, str, str]] = set()
    for coords in coord_map.values():
        for chrom, pos, ref, alt in coords:
            cadd_keys.add((chrom, pos, ref, alt))
    scores = cadd.bulk_lookup(cadd_keys)
    if not scores:
        return

    for a in annotations:
        coords = coord_map.get(a.rsid)
        if not coords:
            continue
        best: float | None = None
        for chrom, pos, ref, alt in coords:
            if a.alt:
                resolved = resolve_strand(a.alt, ref, alt)
                if resolved is None:
                    continue
                score = scores.get((chrom, pos, ref, resolved))
            else:
                score = scores.get((chrom, pos, ref, alt))
            if score is not None and (best is None or score > best):
                best = score
        a.cadd_phred = best


def run_analysis(
    file_path: Path,
    parser: GenotypeParser,
    annotators: list[Annotator],
    skipped_count_provider: Callable[[], int] = lambda: 0,
    *,
    build_override: str | None = None,
    gnomad: GnomadAnnotator | None = None,
    alphamissense: AlphaMissenseAnnotator | None = None,
    cadd: CaddAnnotator | None = None,
) -> AnalysisResult:
    """Stream the file once, query every ready annotator per variant, return results.

    `build_override` short-circuits build detection: when supplied
    (e.g., from `--build grch37`), every variant gets that build and
    the position-data detector is skipped. When None, the pipeline
    buffers the head of the stream until detection is confident, then
    flushes through annotation.

    Annotators are entered into a `contextlib.ExitStack` so their resources
    (e.g., SQLite connections) are deterministically closed.
    """
    metadata = parser.get_metadata(file_path)
    header_build = normalize_build_label(metadata.get("build"))

    annotations: list[Annotation] = []
    total = 0
    diag = _BuildDetectionState(override=build_override, header_build=header_build)

    with contextlib.ExitStack() as stack:
        bound = [stack.enter_context(a) for a in annotators]
        for variant in parser.parse(file_path):
            total += 1
            ready, batch = diag.feed(variant)
            if not ready:
                continue
            for v in batch:
                for annotator in bound:
                    annotations.extend(annotator.annotate(v))
        # End of stream: flush any buffered variants with the best
        # effective build we can resolve (detected → header → default).
        for v in diag.flush():
            for annotator in bound:
                annotations.extend(annotator.annotate(v))

    if gnomad is not None and gnomad.is_ready():
        exact_keys = {(a.rsid, a.alt) for a in annotations if a.alt}
        max_rsids = {a.rsid for a in annotations if not a.alt}
        exact_freq = gnomad.bulk_lookup_by_alt(exact_keys)
        max_freq = gnomad.bulk_lookup(max_rsids)
        for a in annotations:
            if a.alt:
                a.allele_frequency = exact_freq.get((a.rsid, a.alt))
            else:
                a.allele_frequency = max_freq.get(a.rsid)

    if alphamissense is not None and alphamissense.is_ready():
        exact_keys = {(a.rsid, a.alt) for a in annotations if a.alt}
        max_rsids = {a.rsid for a in annotations if not a.alt}
        exact_am = alphamissense.bulk_lookup_by_alt(exact_keys)
        max_am = alphamissense.bulk_lookup(max_rsids)
        for a in annotations:
            hit = exact_am.get((a.rsid, a.alt)) if a.alt else max_am.get(a.rsid)
            if hit is not None:
                a.am_pathogenicity, a.am_class = hit

    if cadd is not None and cadd.is_ready() and gnomad is not None and gnomad.is_ready():
        if getattr(cadd, "_full_mode", False) and diag.effective_build != BUILD_GRCH38:
            logging.getLogger(__name__).warning(
                "CADD full mode requires GRCh38 coordinates; "
                "detected %s — skipping CADD enrichment",
                diag.effective_build,
            )
        else:
            _enrich_cadd(annotations, gnomad, cadd)

    annotators_used = [(a.name, a.version()) for a in annotators]
    if gnomad is not None and gnomad.is_ready():
        annotators_used.append((gnomad.name, gnomad.version()))
    if alphamissense is not None and alphamissense.is_ready():
        annotators_used.append((alphamissense.name, alphamissense.version()))
    if cadd is not None and cadd.is_ready():
        annotators_used.append((cadd.name, cadd.version()))

    return AnalysisResult(
        file_path=file_path,
        parser_name=parser.name,
        parser_display_name=parser.display_name,
        sample_id=metadata["sample_id"],
        build=diag.effective_build,
        total_variants=total,
        skipped_count=skipped_count_provider(),
        annotators_used=annotators_used,
        annotations=annotations,
        build_diagnostics=diag.diagnostics(),
    )


class _BuildDetectionState:
    """Buffer-and-flush state machine for build detection during streaming.

    `feed(variant)` returns (ready, batch). When `ready` is False, the
    variant has been buffered and the caller should keep streaming.
    When True, `batch` contains one or more variants with their build
    field set to the effective build, ready to be annotated.

    `flush()` is called at end of stream to drain anything still
    buffered (which only happens when detection never converged).
    """

    def __init__(self, *, override: str | None, header_build: str | None) -> None:
        self.header_build = header_build
        self.override = override
        # Effective build: starts as override (if given), else None until detection runs.
        self.effective: str | None = override
        self.detected: str | None = None
        self.matched_count = 0
        self.inspected_count = 0
        self._buffer: list[Variant] = []

    @property
    def effective_build(self) -> str:
        """Best-effort effective build at flush time."""
        return self.effective or self.header_build or BUILD_GRCH37

    def feed(self, variant: Variant) -> tuple[bool, list[Variant]]:
        if self.effective is not None:
            return True, [replace(variant, build=self.effective)]
        # Buffering until detection converges or we hit the cap.
        self._buffer.append(variant)
        if variant.rsid in KNOWN_SNP_POSITIONS:
            result = detect_build(self._buffer)
            if result.is_confident:
                self.detected = result.build
                self.matched_count = result.matched
                self.inspected_count = result.inspected
                self.effective = result.build
                batch = [replace(v, build=result.build) for v in self._buffer]
                self._buffer.clear()
                return True, batch
        if len(self._buffer) >= _DETECTION_BUFFER_LIMIT:
            # Buffer full before detection converged. Run partial detection
            # so the GRCh36 safety guard can fire (same logic as flush()).
            result = detect_build(self._buffer)
            if result.build is not None:
                self.detected = result.build
            self.matched_count = result.matched
            self.inspected_count = result.inspected
            if result.build == BUILD_GRCH36:
                self.effective = BUILD_GRCH36
            else:
                self.effective = self.header_build or BUILD_GRCH37
            batch = [replace(v, build=self.effective) for v in self._buffer]
            self._buffer.clear()
            return True, batch
        return False, []

    def flush(self) -> list[Variant]:
        if not self._buffer:
            return []
        # Detection never converged. Re-run on the full buffer to capture
        # partial counts even if not confident.
        result = detect_build(self._buffer)
        if result.is_confident:
            self.detected = result.build
            self.effective = result.build
        else:
            if result.build is not None:
                self.detected = result.build
            # GRCh36 must fail safe: there is no GRCh36 ClinVar cache,
            # so falling back to GRCh37 would silently query wrong
            # coordinates and bypass the GRCh36 safety guard.
            if result.build == BUILD_GRCH36:
                self.effective = BUILD_GRCH36
            else:
                self.effective = self.header_build or BUILD_GRCH37
        self.matched_count = result.matched
        self.inspected_count = result.inspected
        out = [replace(v, build=self.effective) for v in self._buffer]
        self._buffer.clear()
        return out

    def diagnostics(self) -> BuildDiagnostics:
        return BuildDiagnostics(
            header_build=self.header_build,
            detected_build=self.detected,
            effective_build=self.effective_build,
            override=self.override is not None,
            matched_count=self.matched_count,
            inspected_count=self.inspected_count,
        )


__all__ = [
    "BUILD_GRCH37",
    "BUILD_GRCH38",
    "AnalysisResult",
    "BuildDiagnostics",
    "rollup_gwas_duplicates",
    "run_analysis",
]
