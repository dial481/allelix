# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""ClinVar annotator. Source-attributed pathogenicity calls (ADR-0003).

ADR-0021: per-build SQLite caches. ClinVar publishes separate VCFs for
GRCh37 and GRCh38, and the strand orientation of REF/ALT can invert
between builds for the ~0.4% of the genome where the reference
assembly was rebuilt. Carrier-rule matches (ADR-0007) MUST be done
against the same build the user's data is on. The annotator holds one
SQLite cache per build (`clinvar.GRCh37.sqlite`, `clinvar.GRCh38.sqlite`)
and dispatches per-variant by `variant.build`.
"""

from __future__ import annotations

import logging
import sqlite3
from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator, LicenseDescriptor
from allelix.databases import manager as _manager_module
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.manager import (
    download,
    fetch_remote_text,
    get_database_info,
    load_clinvar_vcf,
    stamp_existing_clinvar_cache,
    verify_file_hash,
)
from allelix.models import Annotation

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Variant

logger = logging.getLogger(__name__)

CLINVAR_SUPPORTED_BUILDS: tuple[str, ...] = ("GRCh37", "GRCh38")


def clinvar_db_filename(build: str) -> str:
    """Per-build cache filename. Two coexisting SQLite files per data_dir."""
    return f"clinvar.{build}.sqlite"


def clinvar_record_name(build: str) -> str:
    """`database_versions` row name for a given build."""
    return f"clinvar.{build}"


# Allelix-derived magnitude scoring from ClinVar's CLNSIG. See ADR-0008.
_CLNSIG_MAGNITUDE: dict[str, float] = {
    "pathogenic": 9.0,
    "pathogenic/likely_pathogenic": 8.5,
    "likely_pathogenic": 7.0,
    "drug_response": 6.5,
    "risk_factor": 6.0,
    "uncertain_significance": 4.0,
    "conflicting_interpretations_of_pathogenicity": 4.0,
    "conflicting_classifications_of_pathogenicity": 4.0,
    "not_provided": 2.0,
    "no_classification_for_the_single_variant": 2.0,
    "likely_benign": 2.0,
    "benign/likely_benign": 1.5,
    "benign": 1.0,
}


_BENIGN_CLNSIGS = frozenset({"benign", "likely_benign", "benign/likely_benign"})


def _normalize_clnsig(value: str) -> str:
    return value.strip().lower().replace(" ", "_")


def _magnitude(clnsig: str) -> float:
    return _CLNSIG_MAGNITUDE.get(_normalize_clnsig(clnsig), 5.0)


def _vcf_filename_for_url(url: str) -> str:
    """Pick the right local filename suffix based on the URL."""
    return "clinvar.vcf.gz" if url.endswith(".gz") else "clinvar.vcf"


class ClinVarAnnotator(Annotator):
    """Annotates variants with ClinVar's clinical significance classifications.

    Per-build aware (ADR-0021). At `setup()` time, downloads each
    requested build's VCF (default: both). At `annotate()` time,
    dispatches to the cache matching `variant.build`. If the matching
    cache is missing, the variant is skipped and a warning logged
    (db update needed).
    """

    name: ClassVar[str] = "clinvar"
    display_name: ClassVar[str] = "ClinVar"
    attribution: ClassVar[str] = "ClinVar"
    requires_download: ClassVar[bool] = True
    license: ClassVar[LicenseDescriptor] = LicenseDescriptor(
        spdx="custom-clinvar",
        license_url="https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/",
        attribution_text="ClinVar variant classifications from NCBI.",
        source_url="https://www.ncbi.nlm.nih.gov/clinvar/",
    )

    def __init__(
        self,
        data_dir: Path,
        builds: tuple[str, ...] = CLINVAR_SUPPORTED_BUILDS,
        *,
        include_benign: bool = False,
    ) -> None:
        """Resolve per-build SQLite cache paths within `data_dir`.

        `builds` selects which builds this annotator instance manages.
        Default is both GRCh37 and GRCh38. Passing a single-element
        tuple (e.g. `("GRCh38",)`) restricts setup/refresh to that
        build — used by the CLI's `--build` flag.

        `include_benign` controls whether Benign/Likely_benign annotations
        are emitted. Default False suppresses them (ADR-0008 amendment).
        """
        super().__init__(data_dir)
        self._builds = tuple(builds)
        self._include_benign = include_benign
        for build in self._builds:
            if build not in CLINVAR_SUPPORTED_BUILDS:
                raise ValueError(
                    f"Unsupported ClinVar build {build!r}; expected one of "
                    f"{CLINVAR_SUPPORTED_BUILDS}"
                )
        self._db_paths: dict[str, Path] = {
            build: data_dir / clinvar_db_filename(build) for build in self._builds
        }
        self._conns: dict[str, sqlite3.Connection] = {}
        # ADR-0023: per-build (rsid -> single-base REF) cache. PharmGKB
        # consults this as its primary non-finding filter. Built lazily
        # on first lookup per build.
        self._ref_lookups: dict[str, dict[str, str]] = {}

    def _connection(self, build: str) -> sqlite3.Connection | None:
        """Return a lazy connection to the per-build cache, or None if missing."""
        if build not in self._db_paths:
            return None
        if build not in self._conns:
            db_path = self._db_paths[build]
            if not db_path.exists():
                return None
            self._conns[build] = sqlite3.connect(db_path)
        return self._conns[build]

    def setup(self) -> None:
        """Download each managed build's ClinVar VCF and ingest atomically."""
        for build in self._builds:
            self._setup_one(build)

    def _setup_one(self, build: str) -> None:
        url = _manager_module.CLINVAR_URL_BY_BUILD[build]
        signal = self._fetch_remote_signal_for(build)
        if signal is None:
            msg = (
                f"clinvar ({build}): cannot verify remote freshness signal. "
                "Refresh aborted to avoid persisting an incomplete cache stamp. "
                "Retry, or pass --force if you accept that next `db update` "
                "will re-download to re-establish the signal."
            )
            raise RuntimeError(msg)
        vcf_path = self.data_dir / _vcf_filename_for_url(url)
        download(url, vcf_path)
        try:
            verify_file_hash(vcf_path, "md5", signal.removeprefix("md5:"))
            load_clinvar_vcf(
                vcf_path,
                self._db_paths[build],
                source_url=url,
                remote_signal=signal,
                record_name=clinvar_record_name(build),
            )
        finally:
            try:
                vcf_path.unlink()
            except FileNotFoundError:
                pass
            except OSError:
                logger.warning("Could not remove staged VCF at %s", vcf_path)

    def is_ready(self) -> bool:
        """True iff EVERY managed build has a populated, version-stamped cache.

        Checks ``local_version_tag`` for the current interpreter version.
        Pre-mechanism caches (tag missing or baked into ``remote_signal``)
        are self-healed once via ``stamp_existing_clinvar_cache``.
        """
        for build in self._builds:
            info = get_database_info(self._db_paths[build], clinvar_record_name(build))
            if info is None:
                return False
            tag = info.get("local_version_tag") or ""
            if tag == f"iv:{CLINVAR_INTERPRETER_VERSION}":
                continue
            if stamp_existing_clinvar_cache(self._db_paths[build]):
                continue
            return False
        return True

    def version(self) -> str | None:
        """Composite version string across managed builds.

        Format: `"GRCh37:<v>; GRCh38:<v>"` when both present, or a
        single `<build>:<v>` when only one is managed. None if none.
        """
        parts: list[str] = []
        for build in self._builds:
            info = get_database_info(self._db_paths[build], clinvar_record_name(build))
            if info is not None:
                parts.append(f"{build}:{info['version']}")
        return "; ".join(parts) if parts else None

    def record_count(self) -> int | None:
        """Total record count across managed build caches, or None if none cached."""
        total = 0
        any_present = False
        for build in self._builds:
            info = get_database_info(self._db_paths[build], clinvar_record_name(build))
            if info is not None:
                any_present = True
                total += info["record_count"]
        return total if any_present else None

    def close(self) -> None:
        """Close all open per-build connections. Safe to call repeatedly."""
        for conn in self._conns.values():
            conn.close()
        self._conns.clear()
        self._ref_lookups.clear()

    def reference_for(self, rsid: str, build: str) -> str | None:
        """Return ClinVar's single-base REF allele for `rsid` in `build`, or None.

        ADR-0023: PharmGKB's primary non-finding filter calls this. If the
        return value matches both of the user's alleles, the user is
        homozygous reference and the PharmGKB annotation is a non-finding.

        Lazily builds an in-memory `(rsid -> REF)` map per build on first
        call so subsequent lookups are O(1). Multi-base REFs (indels) are
        skipped — array-based parsers can't call indels, so a multi-base
        REF can't validly suppress a single-base genotype.

        Returns None when ClinVar has no data for the rsid in this build
        (or has only indel REFs). Callers fall through to secondary tiers.
        """
        if build not in self._db_paths:
            return None
        if build not in self._ref_lookups:
            self._ref_lookups[build] = self._load_ref_lookup(build)
        return self._ref_lookups[build].get(rsid)

    def _load_ref_lookup(self, build: str) -> dict[str, str]:
        """Read the per-build cache once and build the `(rsid -> REF)` map."""
        conn = self._connection(build)
        if conn is None:
            return {}
        # Single-base REFs only: indel anchor-base encoding (REF=CTT, etc.)
        # can't suppress a single-base array readout. The per-build cache
        # may have BOTH SNV and indel rows for the same rsid; the WHERE
        # filters those out so we keep only the SNV REF.
        rows = conn.execute(
            "SELECT DISTINCT rsid, ref FROM clinvar_variants WHERE length(ref) = 1"
        ).fetchall()
        out: dict[str, str] = {}
        for rsid, ref in rows:
            # If a rsid has multiple single-base REFs (shouldn't happen at
            # one position but defending against future data shapes), keep
            # the first.
            if rsid not in out:
                out[rsid] = ref
        return out

    def fetch_remote_signal(self) -> str | None:
        r"""Composite freshness signal across managed builds.

        Format: `"GRCh37:md5:<hex>|GRCh38:md5:<hex>"`. Returns None if
        ANY managed build's signal probe fails — the CLI then prints
        "can't verify" and skips refresh per ADR-0012's policy.
        """
        parts: list[str] = []
        for build in self._builds:
            sig = self._fetch_remote_signal_for(build)
            if sig is None:
                return None
            parts.append(f"{build}:{sig}")
        return "|".join(parts) if parts else None

    @staticmethod
    def _fetch_remote_signal_for(build: str) -> str | None:
        body = fetch_remote_text(_manager_module.CLINVAR_URL_BY_BUILD[build] + ".md5")
        if not body:
            return None
        first_token = body.strip().split(None, 1)[0] if body.strip() else ""
        if not first_token:
            return None
        return f"md5:{first_token}"

    def cached_remote_signal(self) -> str | None:
        """Composite cached signal across managed builds. None if any missing."""
        parts: list[str] = []
        for build in self._builds:
            info = get_database_info(self._db_paths[build], clinvar_record_name(build))
            if info is None or info["remote_signal"] is None:
                return None
            sig = info["remote_signal"]
            if not sig:
                return None
            parts.append(f"{build}:{sig}")
        return "|".join(parts) if parts else None

    def annotate(self, variant: Variant) -> list[Annotation]:
        """Return ClinVar annotations whose REF/ALT matches the user's genotype.

        ADR-0007 carrier rule: an entry triggers only if `variant.allele1`
        or `variant.allele2` equals the entry's ALT allele.
        ADR-0011 indel-anchor protection: array-based parsers report
        single-base genotypes; ClinVar's anchor-base indel encoding
        does not match those by string equality.
        ADR-0021: dispatch by `variant.build`. If the matching cache is
        absent, the variant is skipped silently — the user already saw
        the analyze-time build warning.
        """
        if variant.is_no_call:
            return []
        conn = self._connection(variant.build)
        if conn is None:
            return []
        rows = conn.execute(
            "SELECT chromosome, position, ref, alt, clinical_significance, "
            "condition, gene, review_status, allele_id "
            "FROM clinvar_variants WHERE rsid = ?",
            (variant.rsid,),
        ).fetchall()
        annotations: list[Annotation] = []
        carrier_alleles = {variant.allele1, variant.allele2}
        user_is_multibase = len(variant.allele1) > 1 or len(variant.allele2) > 1
        # ADR-0023: report the user's actual diploid call consistently
        # across annotators, not the matched ALT base alone.
        user_diploid = _user_diploid(variant)
        for row in rows:
            (
                _chrom,
                _pos,
                ref,
                alt,
                clnsig,
                condition,
                gene,
                review_status,
                allele_id,
            ) = row
            clinvar_is_indel = len(ref) > 1 or len(alt) > 1
            if clinvar_is_indel and not user_is_multibase:
                continue
            if alt not in carrier_alleles:
                continue
            sig_label = _normalize_clnsig(clnsig) if clnsig else "unknown"
            if not self._include_benign and sig_label in _BENIGN_CLNSIGS:
                continue
            description = (
                f"ClinVar classifies this allele as "
                f"{clnsig.replace('_', ' ') if clnsig else 'unknown significance'}"
            )
            references = [f"clinvar:allele/{allele_id}"] if allele_id else []
            annotations.append(
                Annotation(
                    source=self.name,
                    rsid=variant.rsid,
                    significance=f"clinvar_{sig_label}",
                    category="clinical",
                    magnitude=_magnitude(clnsig),
                    description=description,
                    attribution=self.attribution,
                    genotype_match=user_diploid,
                    references=references,
                    condition=condition or "",
                    gene=gene or "",
                    review_status=review_status or "",
                    alt=alt,
                )
            )
        return annotations


def _user_diploid(variant: Variant) -> str:
    """Render the user's diploid call as a sorted two-letter string.

    Used by ClinVar and PharmGKB so the report's "Genotype" column shows
    the same shape for every annotation regardless of source (ADR-0023).
    SNV: `("G", "A") -> "AG"`. Indel passthrough is verbatim.
    """
    a1, a2 = variant.allele1, variant.allele2
    if len(a1) == 1 and len(a2) == 1:
        return "".join(sorted((a1, a2)))
    return f"{a1}/{a2}"
