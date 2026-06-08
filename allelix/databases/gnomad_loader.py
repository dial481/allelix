# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""gnomAD exome frequency cache loader.

The pre-built SQLite cache is downloaded from HuggingFace during
``db update``. Contains all ~16M exome rsIDs from gnomAD v4.1 with
genomic coordinates (chrom/pos/ref/alt) for future AlphaMissense/CADD
integration.

The cache can also be built locally from gnomAD exome VCFs via
``scripts/build_gnomad_cache.py`` (streaming or local file mode).
"""

from __future__ import annotations

import contextlib
import gzip
import logging
import os
import shutil
import sqlite3
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

logger = logging.getLogger(__name__)

GNOMAD_DB_FILENAME = "gnomad.sqlite"

GNOMAD_CACHE_URL = (
    "https://huggingface.co/datasets/dial481/allelix-gnomad/resolve/main/"
    "exome_frequencies.sqlite.gz"
)


def install_prebuilt_cache(
    gz_path: Path,
    db_path: Path,
    *,
    source_url: str = "",
    remote_signal: str | None = None,
) -> None:
    """Decompress a gzipped pre-built SQLite cache into place.

    The pre-built cache already contains the ``gnomad_frequencies`` table
    and the ``database_versions`` row. This function decompresses and
    stamps the remote signal for freshness tracking.
    """
    tmp_path = db_path.parent / f"{db_path.name}.tmp"
    if tmp_path.exists():
        tmp_path.unlink()

    with gzip.open(gz_path, "rb") as f_in, tmp_path.open("wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    if remote_signal:
        with contextlib.closing(sqlite3.connect(tmp_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET remote_signal = ? WHERE name = 'gnomad'",
                (remote_signal,),
            )
            conn.commit()

    os.replace(tmp_path, db_path)


def cache_exists(db_path: Path) -> bool:
    """Return True if the gnomAD cache exists and has a version row."""
    if not db_path.exists():
        return False
    try:
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            row = conn.execute(
                "SELECT record_count FROM database_versions WHERE name = 'gnomad'"
            ).fetchone()
            return row is not None
    except sqlite3.OperationalError:
        return False
