# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Tests for the config file system."""

from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator, LicenseDescriptor, is_non_commercial
from allelix.config import AllelixConfig, load_config, save_config

if TYPE_CHECKING:
    from pathlib import Path


def _annotator_classes() -> dict[str, type]:
    """Build annotator class map for commercial-mode tests."""
    from allelix.annotators import (
        AlphaMissenseAnnotator,
        ClinVarAnnotator,
        GnomadAnnotator,
        GWASCatalogAnnotator,
        PharmGKBAnnotator,
        SNPediaAnnotator,
    )

    return {
        cls.name: cls
        for cls in [
            ClinVarAnnotator,
            PharmGKBAnnotator,
            GWASCatalogAnnotator,
            SNPediaAnnotator,
            GnomadAnnotator,
            AlphaMissenseAnnotator,
        ]
    }


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
        classes = _annotator_classes()
        assert not config.is_enabled("snpedia", classes)

    def test_commercial_mode_keeps_open_sources(self):
        config = AllelixConfig(commercial=True)
        classes = _annotator_classes()
        assert config.is_enabled("clinvar", classes)
        assert config.is_enabled("pharmgkb", classes)
        assert config.is_enabled("gwas", classes)
        assert config.is_enabled("gnomad", classes)
        assert config.is_enabled("alphamissense", classes)

    def test_commercial_overrides_explicit_enable(self):
        config = AllelixConfig(
            sources={"snpedia": True},
            commercial=True,
        )
        classes = _annotator_classes()
        assert not config.is_enabled("snpedia", classes)

    def test_missing_license_raises_at_definition(self):
        """An annotator subclass without a license ClassVar fails at class definition."""
        import pytest

        with pytest.raises(TypeError, match="must declare a 'license' ClassVar"):

            class BrokenAnnotator(Annotator):
                name: ClassVar[str] = "broken"
                display_name: ClassVar[str] = "Broken"
                attribution: ClassVar[str] = "Broken"

                def setup(self) -> None: ...
                def annotate(self, variant):
                    return []

                def is_ready(self) -> bool:
                    return True

                def version(self) -> str | None:
                    return None

                def close(self) -> None: ...
                def fetch_remote_signal(self) -> str | None:
                    return None

                def cached_remote_signal(self) -> str | None:
                    return None

    def test_commercial_mode_disables_nc_without_explicit_classes(self):
        """Commercial mode blocks NC sources even without passing annotator_classes."""
        config = AllelixConfig(commercial=True)
        assert not config.is_enabled("snpedia")

    def test_nc_spdx_disabled_in_commercial_mode(self):
        """An annotator with a non-commercial SPDX is disabled in commercial mode."""

        class FakeAnnotator(Annotator):
            name: ClassVar[str] = "fake_nc"
            display_name: ClassVar[str] = "Fake NC"
            attribution: ClassVar[str] = "Fake NC"
            license: ClassVar[LicenseDescriptor] = LicenseDescriptor(
                spdx="CC-BY-NC-4.0",
                license_url="https://creativecommons.org/licenses/by-nc/4.0/",
                attribution_text="Fake.",
                source_url="https://example.com",
            )

            def setup(self) -> None: ...
            def annotate(self, variant):
                return []

            def is_ready(self) -> bool:
                return True

            def version(self) -> str | None:
                return None

            def close(self) -> None: ...
            def fetch_remote_signal(self) -> str | None:
                return None

            def cached_remote_signal(self) -> str | None:
                return None

        config = AllelixConfig(commercial=True)
        assert not config.is_enabled("fake_nc", {"fake_nc": FakeAnnotator})
        assert is_non_commercial(FakeAnnotator.license.spdx)


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
        assert not loaded.is_enabled("snpedia", _annotator_classes())

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
