# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Tests for the config file system."""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.config import AllelixConfig, load_config, save_config

if TYPE_CHECKING:
    from pathlib import Path


class TestAllelixConfig:
    def test_default_sources_all_enabled(self):
        config = AllelixConfig()
        assert config.is_enabled("clinvar")
        assert config.is_enabled("pharmgkb")
        assert config.is_enabled("gwas")
        assert config.is_enabled("gnomad")
        assert config.is_enabled("alphamissense")
        assert config.is_enabled("snpedia")

    def test_unknown_source_defaults_enabled(self):
        config = AllelixConfig()
        assert config.is_enabled("future_source")

    def test_disabled_source(self):
        config = AllelixConfig(sources={"gnomad": False})
        assert not config.is_enabled("gnomad")

    def test_commercial_mode_disables_snpedia(self):
        config = AllelixConfig(commercial=True)
        assert not config.is_enabled("snpedia")

    def test_commercial_mode_keeps_open_sources(self):
        config = AllelixConfig(commercial=True)
        assert config.is_enabled("clinvar")
        assert config.is_enabled("pharmgkb")
        assert config.is_enabled("gwas")
        assert config.is_enabled("gnomad")
        assert config.is_enabled("alphamissense")

    def test_commercial_overrides_explicit_enable(self):
        config = AllelixConfig(
            sources={"snpedia": True},
            commercial=True,
        )
        assert not config.is_enabled("snpedia")


class TestLoadSaveConfig:
    def test_creates_default_on_missing(self, tmp_path: Path):
        config = load_config(tmp_path)
        assert config.is_enabled("clinvar")
        assert not config.commercial
        assert (tmp_path / "config.toml").exists()

    def test_roundtrip(self, tmp_path: Path):
        original = AllelixConfig(
            sources={"clinvar": True, "gnomad": False, "snpedia": True},
            commercial=True,
        )
        save_config(tmp_path, original)
        loaded = load_config(tmp_path)
        assert not loaded.is_enabled("gnomad")
        assert loaded.commercial
        assert not loaded.is_enabled("snpedia")

    def test_missing_keys_get_defaults(self, tmp_path: Path):
        (tmp_path / "config.toml").write_text(
            "[sources]\nclinvar = false\n\n[license]\ncommercial = false\n"
        )
        config = load_config(tmp_path)
        assert not config.is_enabled("clinvar")
        assert config.is_enabled("gnomad")

    def test_corrupt_values_ignored(self, tmp_path: Path):
        (tmp_path / "config.toml").write_text(
            '[sources]\nclinvar = "yes"\n\n[license]\ncommercial = false\n'
        )
        config = load_config(tmp_path)
        assert config.is_enabled("clinvar")

    def test_serialized_format(self, tmp_path: Path):
        config = AllelixConfig(commercial=True)
        save_config(tmp_path, config)
        text = (tmp_path / "config.toml").read_text()
        assert "[sources]" in text
        assert "[license]" in text
        assert "commercial = true" in text
