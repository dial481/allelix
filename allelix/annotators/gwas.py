# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""GWAS Catalog annotator. Source-attributed trait associations."""

from __future__ import annotations

import logging
import sqlite3
import zipfile
from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator
from allelix.databases.gwas_loader import (
    GWAS_CATALOG_URL,
    GWAS_DB_FILENAME,
    load_gwas_tsv,
    schema_is_current,
)
from allelix.databases.manager import (
    download,
    get_database_info,
    head_request_headers,
)
from allelix.models import Annotation

_EXCLUDED_TRAIT_CATEGORIES = frozenset(
    {
        "body_measurement",
        "lipid_measurement",
        "hematological_measurement",
        "other_measurement",
        "behavioral",
    }
)

_MUST_INCLUDE_RSIDS = frozenset(
    {
        "rs10737680",  # CFH — age-related macular degeneration
        "rs11209026",  # IL23R — inflammatory bowel disease
        "rs9271366",  # HLA-DRB1 — multiple sclerosis
    }
)

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Variant

logger = logging.getLogger(__name__)


def _magnitude(p_value: float | None, or_beta: float | None) -> float:
    """Derive magnitude from p-value and optional effect size."""
    if p_value is None:
        base = 2.0
    elif p_value < 5e-100:
        base = 8.0
    elif p_value < 5e-20:
        base = 7.0
    elif p_value < 5e-8:
        base = 6.0
    elif p_value < 5e-6:
        base = 4.0
    elif p_value < 5e-4:
        base = 3.0
    else:
        base = 2.0

    if or_beta is not None and or_beta > 0:
        if or_beta >= 3.0 or or_beta <= 0.33:
            base = min(base + 1.0, 9.0)
        elif or_beta >= 2.0 or or_beta <= 0.5:
            base = min(base + 0.5, 9.0)

    return base


_UNKNOWN_RISK_ALLELE_MAG_CAP = 3.0


class GWASCatalogAnnotator(Annotator):
    """Annotates variants with GWAS Catalog trait associations."""

    name: ClassVar[str] = "gwas"
    display_name: ClassVar[str] = "GWAS Catalog"
    attribution: ClassVar[str] = "GWAS Catalog"
    requires_download: ClassVar[bool] = True

    def __init__(self, data_dir: Path, *, filter_traits: bool = True) -> None:
        """Initialize with path to the data directory."""
        super().__init__(data_dir)
        self._db_path = data_dir / GWAS_DB_FILENAME
        self._conn: sqlite3.Connection | None = None
        self._filter_traits = filter_traits

    def _connection(self) -> sqlite3.Connection:
        if self._conn is None:
            self._conn = sqlite3.connect(self._db_path)
        return self._conn

    def setup(self) -> None:
        """Download GWAS Catalog associations ZIP, extract TSV, and ingest."""
        url = GWAS_CATALOG_URL
        signal = self.fetch_remote_signal()
        if signal is None:
            msg = (
                "gwas: cannot verify remote freshness signal. "
                "Refresh aborted to avoid persisting an incomplete cache stamp. "
                "Retry, or pass --force if you accept that next `db update` "
                "will re-download to re-establish the signal."
            )
            raise RuntimeError(msg)
        zip_path = self.data_dir / "gwas_catalog_associations.zip"
        tsv_path = self.data_dir / "gwas_catalog_associations.tsv"
        download(url, zip_path)
        try:
            with zipfile.ZipFile(zip_path) as zf:
                tsv_names = [n for n in zf.namelist() if n.endswith(".tsv")]
                if not tsv_names:
                    msg = f"No .tsv file found in {zip_path}"
                    raise RuntimeError(msg)
                zf.extract(tsv_names[0], self.data_dir)
                extracted = self.data_dir / tsv_names[0]
                if extracted != tsv_path:
                    extracted.rename(tsv_path)
            load_gwas_tsv(tsv_path, self._db_path, source_url=url, remote_signal=signal)
        finally:
            try:
                zip_path.unlink()
            except FileNotFoundError:
                pass
            except OSError:
                logger.warning("Could not remove staged file at %s", zip_path)

    def is_ready(self) -> bool:
        """Return True when the local GWAS cache exists and has current schema.

        When the categorizer version has bumped and the raw TSV is still
        on disk (retained since the last ``db update``), auto-reingests
        from the cached TSV without re-downloading.
        """
        if get_database_info(self._db_path, "gwas") is None:
            return False
        if schema_is_current(self._db_path):
            return True
        tsv_path = self.data_dir / "gwas_catalog_associations.tsv"
        if tsv_path.exists():
            print(
                "GWAS categorizer changed — re-ingesting from cached TSV...",
                flush=True,
            )
            try:
                load_gwas_tsv(
                    tsv_path,
                    self._db_path,
                    source_url=GWAS_CATALOG_URL,
                    remote_signal=self.cached_remote_signal(),
                )
            except Exception:
                logger.warning("Auto-reingest from cached TSV failed", exc_info=True)
                return False
            return schema_is_current(self._db_path)
        return False

    def version(self) -> str | None:
        """Return the cached database version string, or None."""
        info = get_database_info(self._db_path, "gwas")
        return info["version"] if info else None

    def record_count(self) -> int | None:
        """Return the number of cached GWAS association records, or None."""
        info = get_database_info(self._db_path, "gwas")
        return info["record_count"] if info else None

    def close(self) -> None:
        """Close the SQLite connection if open."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def fetch_remote_signal(self) -> str | None:
        """Probe the GWAS Catalog URL for ETag or Last-Modified."""
        headers = head_request_headers(GWAS_CATALOG_URL)
        if headers is None:
            return None
        etag = headers.get("ETag") or headers.get("Etag")
        last_modified = headers.get("Last-Modified") or headers.get("Last-modified")
        if etag:
            return f"etag:{etag.strip()}"
        if last_modified:
            return f"lm:{last_modified.strip()}"
        return None

    def cached_remote_signal(self) -> str | None:
        """Return the remote signal stored during the last successful ingest."""
        info = get_database_info(self._db_path, "gwas")
        if not info or not info["remote_signal"]:
            return None
        sig = info["remote_signal"]
        return sig.split("|cv:")[0] or None

    def annotate(self, variant: Variant) -> list[Annotation]:
        """Return GWAS Catalog annotations for variants the user carries.

        Carrier matching uses the risk allele when specified. When the
        GWAS entry doesn't specify a risk allele, the annotation fires
        on rsid match alone with a magnitude penalty.
        """
        if variant.is_no_call:
            return []

        sql = (
            "SELECT risk_allele, trait, p_value, or_beta, gene, "
            "study_accession, pubmed_id, trait_category "
            "FROM gwas_associations WHERE rsid = ?"
        )
        rows = self._connection().execute(sql, (variant.rsid,)).fetchall()
        annotations: list[Annotation] = []
        user_diploid = _user_diploid(variant)

        for row in rows:
            (
                risk_allele,
                trait,
                p_value,
                or_beta,
                gene,
                study_accession,
                pubmed_id,
                trait_category,
            ) = row

            if self._filter_traits and trait_category in _EXCLUDED_TRAIT_CATEGORIES:
                continue

            if risk_allele is not None:
                if variant.allele1 != risk_allele and variant.allele2 != risk_allele:
                    continue
                mag = _magnitude(p_value, or_beta)
                risk_note = ""
            else:
                # ADR-0024: unknown risk allele fires on rsID match alone
                # but capped at 3.0 so it doesn't pass typical --min-magnitude
                # thresholds. Without knowing which allele is the risk allele,
                # we can't apply the carrier rule (ADR-0007).
                mag = min(_magnitude(p_value, or_beta), _UNKNOWN_RISK_ALLELE_MAG_CAP)
                risk_note = " (risk allele not specified in study)"

            p_str = f"p={p_value:.1e}" if p_value is not None else "p=N/A"
            gene_str = gene or "—"
            description = f"GWAS Catalog: {trait} ({p_str}, gene: {gene_str}){risk_note}"

            references: list[str] = []
            if pubmed_id:
                references.append(f"pubmed:{pubmed_id}")
            if study_accession:
                references.append(f"gwas:{study_accession}")

            annotations.append(
                Annotation(
                    source=self.name,
                    rsid=variant.rsid,
                    significance="gwas_association",
                    category="trait",
                    magnitude=mag,
                    description=description,
                    attribution=self.attribution,
                    genotype_match=user_diploid,
                    references=references,
                    condition=trait,
                    gene=gene or "",
                    is_must_include=variant.rsid in _MUST_INCLUDE_RSIDS,
                )
            )
        return annotations


def _user_diploid(variant: Variant) -> str:
    """Sorted two-letter diploid for SNVs; indel passthrough verbatim."""
    a1, a2 = variant.allele1, variant.allele2
    if len(a1) == 1 and len(a2) == 1:
        return "".join(sorted((a1, a2)))
    return f"{a1}/{a2}"
