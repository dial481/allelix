# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""AlphaMissense variant pathogenicity enrichment.

AlphaMissense is not a clinical annotator — it does not produce
Annotation objects. It enriches existing annotations with missense
variant pathogenicity predictions. The pipeline calls
``bulk_lookup()`` after all annotators have run, and stamps each
annotation's ``am_pathogenicity`` and ``am_class`` fields.

License: CC BY 4.0. Attribution: Cheng et al., Science 2023
(doi:10.1126/science.adg7492).
"""

from __future__ import annotations

import logging
import sqlite3
from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator
from allelix.databases.alphamissense_loader import (
    ALPHAMISSENSE_CACHE_URL,
    ALPHAMISSENSE_DB_FILENAME,
    install_prebuilt_cache,
)
from allelix.databases.manager import download, get_database_info, head_request_headers

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Annotation, Variant

logger = logging.getLogger(__name__)

_BULK_BATCH_SIZE = 900


class AlphaMissenseAnnotator(Annotator):
    """Missense variant pathogenicity enrichment from AlphaMissense.

    Subclasses Annotator for ``db update`` / ``db status`` / ``is_ready()``
    integration. ``annotate()`` always returns ``[]`` — AlphaMissense does
    not participate in the per-variant annotation loop.
    """

    name: ClassVar[str] = "alphamissense"
    display_name: ClassVar[str] = "AlphaMissense"
    attribution: ClassVar[str] = "AlphaMissense"
    requires_download: ClassVar[bool] = True

    def __init__(self, data_dir: Path) -> None:
        """Bind to the data directory."""
        super().__init__(data_dir)
        self._db_path = data_dir / ALPHAMISSENSE_DB_FILENAME
        self._conn: sqlite3.Connection | None = None

    def _connection(self) -> sqlite3.Connection:
        if self._conn is None:
            if not self._db_path.exists():
                raise FileNotFoundError(
                    f"AlphaMissense cache not found at {self._db_path}. "
                    "Run `allelix db update` first."
                )
            self._conn = sqlite3.connect(self._db_path)
        return self._conn

    def setup(self) -> None:
        """Download the pre-built AlphaMissense cache from HuggingFace."""
        signal = self.fetch_remote_signal()
        gz_path = self.data_dir / "alphamissense.sqlite.gz"
        download(ALPHAMISSENSE_CACHE_URL, gz_path)
        install_prebuilt_cache(
            gz_path,
            self._db_path,
            source_url=ALPHAMISSENSE_CACHE_URL,
            remote_signal=signal,
        )
        try:
            gz_path.unlink()
        except OSError:
            logger.warning("Could not remove staged file at %s", gz_path)

    def is_ready(self) -> bool:
        """True when the AlphaMissense SQLite cache exists and is queryable."""
        return get_database_info(self._db_path, "alphamissense") is not None

    def version(self) -> str | None:
        """Return the cached database version, or None."""
        info = get_database_info(self._db_path, "alphamissense")
        return info["version"] if info else None

    def record_count(self) -> int | None:
        """Return the number of variants in the cache, or None."""
        info = get_database_info(self._db_path, "alphamissense")
        return info["record_count"] if info else None

    def close(self) -> None:
        """Close the SQLite connection if open."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def fetch_remote_signal(self) -> str | None:
        """Probe the HuggingFace asset URL for ETag or Last-Modified."""
        headers = head_request_headers(ALPHAMISSENSE_CACHE_URL)
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
        info = get_database_info(self._db_path, "alphamissense")
        if not info or not info["remote_signal"]:
            return None
        return info["remote_signal"]

    def annotate(self, variant: Variant) -> list[Annotation]:
        """Not used — AlphaMissense enriches, does not annotate. Always returns []."""
        return []

    def lookup(self, rsid: str) -> tuple[float, str] | None:
        """Return (am_pathogenicity, am_class) for a single rsID, or None."""
        conn = self._connection()
        row = conn.execute(
            "SELECT MAX(am_pathogenicity), am_class FROM alphamissense_scores WHERE rsid = ?",
            (rsid,),
        ).fetchone()
        if row is None or row[0] is None:
            return None
        return (row[0], row[1])

    def bulk_lookup(self, rsids: set[str]) -> dict[str, tuple[float, str]]:
        """Return ``{rsid: (am_pathogenicity, am_class)}`` for found rsIDs.

        When an rsID maps to multiple rows (different ref/alt at the same
        position), returns the highest pathogenicity score and its class —
        consistent with gnomAD's ``MAX(af)`` approach.

        Batches into chunks of 900 to stay within SQLite's variable limit.
        """
        if not rsids:
            return {}
        conn = self._connection()
        result: dict[str, tuple[float, str]] = {}
        rsid_list = list(rsids)
        for i in range(0, len(rsid_list), _BULK_BATCH_SIZE):
            batch = rsid_list[i : i + _BULK_BATCH_SIZE]
            placeholders = ",".join("?" * len(batch))
            rows = conn.execute(
                f"SELECT rsid, MAX(am_pathogenicity), am_class"
                f" FROM alphamissense_scores"
                f" WHERE rsid IN ({placeholders}) GROUP BY rsid",
                batch,
            ).fetchall()
            for rsid, score, cls in rows:
                if score is not None:
                    result[rsid] = (score, cls)
        return result
