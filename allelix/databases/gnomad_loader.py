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
    "https://huggingface.co/datasets/dial481/allelix-gnomad/resolve/main/gnomad.sqlite.gz"
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
    gz_size = gz_path.stat().st_size
    free = shutil.disk_usage(db_path.parent).free
    needed = gz_size * 6
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

    from allelix.databases.manager import _ensure_local_version_tag_column

    with contextlib.closing(sqlite3.connect(tmp_path)) as conn:
        if remote_signal:
            from allelix.databases.manager import stamp_remote_signal

            stamp_remote_signal(conn, "gnomad", remote_signal, source_url)
        from allelix.annotators._versions import GNOMAD_SCHEMA_VERSION

        _ensure_local_version_tag_column(conn)
        conn.execute(
            "UPDATE database_versions SET local_version_tag = ? WHERE name = 'gnomad'",
            (f"sv:{GNOMAD_SCHEMA_VERSION}",),
        )
        conn.commit()

    os.replace(tmp_path, db_path)
