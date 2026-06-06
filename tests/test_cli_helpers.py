# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Unit tests for private CLI helpers."""

from __future__ import annotations

from allelix.cli import _chrom_sort_key, _percent


class TestChromSortKey:
    def test_autosomes_then_sex_then_mt(self):
        keys = sorted(["MT", "X", "1", "22", "Y", "2"], key=_chrom_sort_key)
        assert keys == ["1", "2", "22", "X", "Y", "MT"]

    def test_unknown_chrom_falls_to_alphabetical(self):
        keys = sorted(["X", "1", "WEIRD", "2"], key=_chrom_sort_key)
        assert keys == ["1", "2", "X", "WEIRD"]


class TestPercent:
    def test_zero_total_does_not_divide(self):
        assert _percent(0, 0) == "0.00%"

    def test_basic(self):
        assert _percent(1, 2) == "50.00%"
