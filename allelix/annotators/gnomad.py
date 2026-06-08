# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""gnomAD population frequency enrichment.

gnomAD is not a clinical annotator — it does not produce Annotation
objects. It enriches existing annotations with population allele
frequency context. The pipeline calls ``bulk_lookup()`` after all
annotators have run, and stamps each annotation's ``allele_frequency``
field.

License: ODbL v1.0 (Open Database License). We extract only rsID +
allele frequencies (no SpliceAI or other restrictively licensed fields).
"""

from __future__ import annotations

import logging
import sqlite3
from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator
from allelix.databases.gnomad_loader import (
    GNOMAD_CACHE_URL,
    GNOMAD_DB_FILENAME,
    install_prebuilt_cache,
)
from allelix.databases.manager import download, get_database_info, head_request_headers

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Annotation, Variant

logger = logging.getLogger(__name__)

_BULK_BATCH_SIZE = 900


class GnomadAnnotator(Annotator):
    """Population frequency enrichment from gnomAD.

    Subclasses Annotator for ``db update`` / ``db status`` / ``is_ready()``
    integration. ``annotate()`` always returns ``[]`` — gnomAD does not
    participate in the per-variant annotation loop.
    """

    name: ClassVar[str] = "gnomad"
    display_name: ClassVar[str] = "gnomAD"
    attribution: ClassVar[str] = "gnomAD"
    requires_download: ClassVar[bool] = True

    def __init__(self, data_dir: Path) -> None:
        """Bind to the data directory."""
        super().__init__(data_dir)
        self._db_path = data_dir / GNOMAD_DB_FILENAME
        self._conn: sqlite3.Connection | None = None

    def _connection(self) -> sqlite3.Connection:
        if self._conn is None:
            self._conn = sqlite3.connect(self._db_path)
        return self._conn

    def setup(self) -> None:
        """Download the pre-built gnomAD exome frequency cache from HuggingFace."""
        signal = self.fetch_remote_signal()
        gz_path = self.data_dir / "gnomad_exome_frequencies.sqlite.gz"
        download(GNOMAD_CACHE_URL, gz_path)
        install_prebuilt_cache(
            gz_path,
            self._db_path,
            source_url=GNOMAD_CACHE_URL,
            remote_signal=signal,
        )
        try:
            gz_path.unlink()
        except OSError:
            logger.warning("Could not remove staged file at %s", gz_path)

    def is_ready(self) -> bool:
        """True when the gnomAD SQLite cache exists and is queryable."""
        return get_database_info(self._db_path, "gnomad") is not None

    def version(self) -> str | None:
        """Return the cached database version, or None."""
        info = get_database_info(self._db_path, "gnomad")
        return info["version"] if info else None

    def record_count(self) -> int | None:
        """Return the number of rsIDs in the cache, or None."""
        info = get_database_info(self._db_path, "gnomad")
        return info["record_count"] if info else None

    def close(self) -> None:
        """Close the SQLite connection if open."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def fetch_remote_signal(self) -> str | None:
        """Probe the HuggingFace asset URL for ETag or Last-Modified."""
        headers = head_request_headers(GNOMAD_CACHE_URL)
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
        """Return the remote signal stored at last successful download."""
        info = get_database_info(self._db_path, "gnomad")
        if not info or not info["remote_signal"]:
            return None
        return info["remote_signal"]

    def annotate(self, variant: Variant) -> list[Annotation]:
        """Not used — gnomAD enriches, does not annotate. Always returns []."""
        return []

    def lookup(self, rsid: str) -> float | None:
        """Return global allele frequency for a single rsID, or None."""
        conn = self._connection()
        row = conn.execute(
            "SELECT MAX(af) FROM gnomad_frequencies WHERE rsid = ?", (rsid,)
        ).fetchone()
        return row[0] if row else None

    def bulk_lookup(self, rsids: set[str]) -> dict[str, float]:
        """Return ``{rsid: af}`` for all rsIDs found in the cache.

        Batches into chunks of 900 to stay within SQLite's variable limit.
        """
        if not rsids:
            return {}
        conn = self._connection()
        result: dict[str, float] = {}
        rsid_list = list(rsids)
        for i in range(0, len(rsid_list), _BULK_BATCH_SIZE):
            batch = rsid_list[i : i + _BULK_BATCH_SIZE]
            placeholders = ",".join("?" * len(batch))
            rows = conn.execute(
                f"SELECT rsid, MAX(af) FROM gnomad_frequencies"
                f" WHERE rsid IN ({placeholders}) GROUP BY rsid",
                batch,
            ).fetchall()
            for rsid, af in rows:
                if af is not None:
                    result[rsid] = af
        return result
