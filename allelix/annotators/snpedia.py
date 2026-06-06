# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""SNPedia annotator. Structured SQL lookups against pre-parsed genotype data.

Reads from the ``snpedia_genotypes`` table in the SNPedia SQLite archive.
Raw wiki markup is scraped by ``scripts/scrape_snpedia.py`` into ``pages``.
On first use, this annotator automatically parses raw markup into structured
columns (one-time operation). After that, all queries are pure SQL.

SNPedia content is CC-BY-NC-SA 3.0 US. Attribution is required in all
reports.
"""

from __future__ import annotations

import contextlib
import logging
import sqlite3
from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator, is_clinvar_homref
from allelix.models import Annotation

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path

    from allelix.models import Variant

logger = logging.getLogger(__name__)

SNPEDIA_DB_FILENAME = "snpedia.sqlite"
SNPEDIA_RECORD_NAME = "snpedia"

_REPUTE_CATEGORY: dict[str, str] = {
    "good": "trait",
    "bad": "clinical",
    "not set": "trait",
    "": "trait",
}

_SUMMARY_SUPPRESS_SUBSTRINGS: tuple[str, ...] = (
    "mis-oriented",
    "mis-orientation",
    "wrong strand",
    "orientation uncertain",
)


class SNPediaAnnotator(Annotator):
    """Annotates variants with SNPedia genotype data via structured SQL lookups."""

    name: ClassVar[str] = "snpedia"
    display_name: ClassVar[str] = "SNPedia"
    attribution: ClassVar[str] = "SNPedia"
    requires_download: ClassVar[bool] = False

    def __init__(
        self,
        data_dir: Path,
        clinvar_ref_provider: Callable[[str, str], str | None] | None = None,
    ) -> None:
        """Initialize with path to the data directory.

        ``clinvar_ref_provider`` is a ``(rsid, build) -> ref_base | None``
        callable used by the ADR-0023 hom-ref check. In production it is
        wired to ``ClinVarAnnotator.reference_for``. ``None`` disables the
        check (tests, standalone use).
        """
        super().__init__(data_dir)
        self._db_path = data_dir / SNPEDIA_DB_FILENAME
        self._conn: sqlite3.Connection | None = None
        self._clinvar_ref_provider = clinvar_ref_provider

    def _connection(self) -> sqlite3.Connection:
        if self._conn is None:
            self._conn = sqlite3.connect(self._db_path)
        return self._conn

    def setup(self) -> None:
        """No-op. The SNPedia database is populated by scripts/scrape_snpedia.py."""

    def is_ready(self) -> bool:
        """Return True when the parsed SNPedia genotype table exists and has data.

        If raw pages exist but the structured table does not, automatically
        parses the raw markup (one-time operation, ~2 minutes).
        """
        if not self._db_path.exists():
            return False
        try:
            from allelix.databases.snpedia_parser import (
                detect_raw_table,
                parse_raw_pages,
                parser_is_current,
            )

            with contextlib.closing(sqlite3.connect(self._db_path)) as conn:
                has_rows = False
                with contextlib.suppress(sqlite3.OperationalError):
                    has_rows = (
                        conn.execute("SELECT COUNT(*) FROM snpedia_genotypes").fetchone()[0] > 0
                    )

                needs_reparse = has_rows and not parser_is_current(conn)
                if has_rows and not needs_reparse:
                    return True

                raw_table = detect_raw_table(conn)
                if raw_table is None:
                    return False

                snp_count = conn.execute(
                    f"SELECT COUNT(*) FROM {raw_table} WHERE category='snp'"
                ).fetchone()[0]
                genotype_count = conn.execute(
                    f"SELECT COUNT(*) FROM {raw_table} WHERE category='genotype'"
                ).fetchone()[0]

            reason = "parser version changed" if needs_reparse else "one-time"
            print(
                f"Parsing {snp_count} SNP pages + {genotype_count} genotype pages"
                f" into structured table ({reason}, ~5 min)...",
                flush=True,
            )
            parsed = parse_raw_pages(str(self._db_path))
            print(f"Parsed {parsed} SNPedia genotype rows.", flush=True)
            return parsed > 0
        except (sqlite3.OperationalError, sqlite3.DatabaseError):
            return False

    def version(self) -> str | None:
        """Return a version string from the database_versions table."""
        if not self._db_path.exists():
            return None
        try:
            with contextlib.closing(sqlite3.connect(self._db_path)) as conn:
                row = conn.execute(
                    "SELECT version FROM database_versions WHERE name = ?",
                    (SNPEDIA_RECORD_NAME,),
                ).fetchone()
                if row and row[0]:
                    return row[0]
                return None
        except (sqlite3.OperationalError, sqlite3.DatabaseError):
            return None

    def record_count(self) -> int | None:
        """Return the number of genotype rows in the structured table."""
        if not self._db_path.exists():
            return None
        try:
            with contextlib.closing(sqlite3.connect(self._db_path)) as conn:
                count = conn.execute("SELECT COUNT(*) FROM snpedia_genotypes").fetchone()[0]
                return count
        except (sqlite3.OperationalError, sqlite3.DatabaseError):
            return None

    def close(self) -> None:
        """Close the SQLite connection."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def fetch_remote_signal(self) -> str | None:
        """SNPedia is frozen — no remote freshness signal."""
        return None

    def cached_remote_signal(self) -> str | None:
        """No remote signal for a locally-scraped archive."""
        return None

    def annotate(self, variant: Variant) -> list[Annotation]:
        """Return SNPedia annotations matching the user's genotype."""
        if variant.is_no_call:
            return []

        snp_id = variant.rsid.lower()
        if snp_id.startswith("rs"):
            snp_num = snp_id[2:]
            snp_url_path = f"Rs{snp_num}"
        elif snp_id.startswith("i"):
            snp_num = snp_id[1:]
            snp_url_path = f"I{snp_num}"
        else:
            return []

        if not snp_num or not snp_num.isdigit():
            return []

        if snp_id.startswith("rs") and is_clinvar_homref(variant, self._clinvar_ref_provider):
            return []

        a1, a2 = variant.allele1.upper(), variant.allele2.upper()
        sorted_alleles = (a1, a2) if a1 <= a2 else (a2, a1)

        conn = self._connection()
        rows = conn.execute(
            "SELECT allele1, allele2, magnitude, repute, summary, gene "
            "FROM snpedia_genotypes "
            "WHERE rsid = ? AND allele1 = ? AND allele2 = ?",
            (snp_id, sorted_alleles[0], sorted_alleles[1]),
        ).fetchall()

        annotations: list[Annotation] = []
        for allele1, allele2, magnitude, repute, summary, gene in rows:
            if not summary:
                continue

            summary_lower = summary.lower()
            if any(p in summary_lower for p in _SUMMARY_SUPPRESS_SUBSTRINGS):
                continue

            if magnitude is None:
                magnitude = 0.0

            repute_lower = (repute or "").strip().lower()
            category = _REPUTE_CATEGORY.get(repute_lower, "trait")

            description = f"SNPedia: {summary}"
            genotype_match = f"{allele1}{allele2}"

            annotations.append(
                Annotation(
                    source=self.name,
                    rsid=variant.rsid,
                    significance=f"snpedia_{repute_lower}" if repute_lower else "snpedia_genotype",
                    category=category,
                    magnitude=magnitude,
                    description=description,
                    attribution=self.attribution,
                    genotype_match=genotype_match,
                    references=[f"https://www.snpedia.com/index.php/{snp_url_path}"],
                    gene=gene or "",
                )
            )

        return annotations
