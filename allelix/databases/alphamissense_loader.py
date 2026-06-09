# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""AlphaMissense pathogenicity cache loader.

The pre-built SQLite cache is downloaded from HuggingFace during
``db update``. Contains 71M missense variant scores from AlphaMissense
with genomic coordinates and rsIDs (joined from gnomAD at build time).

The cache can also be built locally from the Zenodo TSV via
``scripts/build_alphamissense_cache.py``.
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

ALPHAMISSENSE_DB_FILENAME = "alphamissense.sqlite"

ALPHAMISSENSE_CACHE_URL = (
    "https://huggingface.co/datasets/dial481/allelix-alphamissense/resolve/main/"
    "alphamissense.sqlite.gz"
)


def install_prebuilt_cache(
    gz_path: Path,
    db_path: Path,
    *,
    source_url: str = "",
    remote_signal: str | None = None,
) -> None:
    """Decompress a gzipped pre-built SQLite cache into place."""
    gz_size = gz_path.stat().st_size
    free = shutil.disk_usage(db_path.parent).free
    needed = gz_size * 5
    if free < needed:
        raise OSError(
            f"Not enough disk space to decompress {gz_path.name}: "
            f"{free / 1e9:.1f} GB free, need ~{needed / 1e9:.1f} GB. "
            "Free up space and retry."
        )

    tmp_path = db_path.parent / f"{db_path.name}.tmp"
    if tmp_path.exists():
        tmp_path.unlink()

    with gzip.open(gz_path, "rb") as f_in, tmp_path.open("wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    if remote_signal:
        with contextlib.closing(sqlite3.connect(tmp_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET remote_signal = ? WHERE name = 'alphamissense'",
                (remote_signal,),
            )
            conn.commit()

    os.replace(tmp_path, db_path)
