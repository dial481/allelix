# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""User configuration for source management and license mode."""

from __future__ import annotations

import logging
import tomllib
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

logger = logging.getLogger(__name__)

CONFIG_FILENAME = "config.toml"

NON_COMMERCIAL_SOURCES: frozenset[str] = frozenset({"snpedia"})

_DEFAULT_SOURCES: dict[str, bool] = {
    "clinvar": True,
    "pharmgkb": True,
    "gwas": True,
    "gnomad": True,
    "alphamissense": True,
    "snpedia": True,
}


@dataclass
class AllelixConfig:
    """Persistent user configuration loaded from ``config.toml``."""

    sources: dict[str, bool] = field(default_factory=lambda: dict(_DEFAULT_SOURCES))
    commercial: bool = False

    def is_enabled(self, source_name: str) -> bool:
        """Check if a source is enabled, respecting commercial mode."""
        if self.commercial and source_name in NON_COMMERCIAL_SOURCES:
            return False
        return self.sources.get(source_name, True)


def _config_path(data_dir: Path) -> Path:
    return data_dir / CONFIG_FILENAME


def _serialize(config: AllelixConfig) -> str:
    """Serialize config to TOML string. Manual formatting — no dependency."""
    lines = ["[sources]"]
    for name, enabled in sorted(config.sources.items()):
        lines.append(f"{name} = {str(enabled).lower()}")
    lines.append("")
    lines.append("[license]")
    lines.append(f"commercial = {str(config.commercial).lower()}")
    lines.append("")
    return "\n".join(lines)


def save_config(data_dir: Path, config: AllelixConfig) -> None:
    """Write config to disk."""
    path = _config_path(data_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_serialize(config), encoding="utf-8")


def load_config(data_dir: Path) -> AllelixConfig:
    """Load config from disk. Creates default config if absent."""
    path = _config_path(data_dir)
    if not path.exists():
        config = AllelixConfig()
        save_config(data_dir, config)
        logger.info("Created default config at %s", path)
        return config

    with open(path, "rb") as fh:
        raw = tomllib.load(fh)

    sources = dict(_DEFAULT_SOURCES)
    if "sources" in raw and isinstance(raw["sources"], dict):
        for key, val in raw["sources"].items():
            if isinstance(val, bool):
                sources[key] = val

    commercial = False
    if "license" in raw and isinstance(raw["license"], dict):
        commercial = bool(raw["license"].get("commercial", False))

    return AllelixConfig(sources=sources, commercial=commercial)
