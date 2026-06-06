# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Tests for strand-flip / complement / ambiguity helpers."""

from __future__ import annotations

from allelix.utils.allele import complement, flip_genotype, is_strand_ambiguous


class TestComplement:
    def test_single_bases(self):
        assert complement("A") == "T"
        assert complement("T") == "A"
        assert complement("C") == "G"
        assert complement("G") == "C"

    def test_no_call_unchanged(self):
        assert complement("-") == "-"
        assert complement("") == ""

    def test_unknown_letter_unchanged(self):
        assert complement("N") == "N"

    def test_multibase_indel_reverses_and_complements(self):
        # CTT (forward) → AAG (reverse complement)
        assert complement("CTT") == "AAG"
        assert complement("AAG") == "CTT"


class TestFlipGenotype:
    def test_diploid_flip(self):
        assert flip_genotype("C", "T") == ("G", "A")

    def test_no_call_preserved(self):
        assert flip_genotype("-", "A") == ("-", "T")


class TestIsStrandAmbiguous:
    def test_at_pair_is_ambiguous(self):
        assert is_strand_ambiguous("A", "T")
        assert is_strand_ambiguous("T", "A")

    def test_cg_pair_is_ambiguous(self):
        assert is_strand_ambiguous("C", "G")
        assert is_strand_ambiguous("G", "C")

    def test_normal_pair_not_ambiguous(self):
        assert not is_strand_ambiguous("A", "G")
        assert not is_strand_ambiguous("C", "T")

    def test_indel_not_ambiguous(self):
        assert not is_strand_ambiguous("CTT", "C")

    def test_unknown_letter_not_ambiguous(self):
        assert not is_strand_ambiguous("A", "N")
