# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Shared helpers for parsers with CSV or concatenated-genotype formats.

Used by FTDNA, MyHeritage, and Living DNA parsers. Extracted here to avoid
duplicating the genotype-splitting and CSV-line-splitting logic across
structurally similar formats.
"""

from __future__ import annotations

import logging

from allelix.models import NO_CALL_MARKER

logger = logging.getLogger(__name__)


def split_csv_line(line: str) -> list[str]:
    """Split a comma-delimited line and strip double-quotes from each field.

    Handles single-quoted, double-quoted, and double-double-quoted fields
    (the MyHeritage "extra quotes" variant).
    """
    return [field.strip().strip('"') for field in line.split(",")]


def split_genotype(genotype: str) -> tuple[str, str]:
    """Split a concatenated genotype field into two alleles.

    ``"AG"`` -> ``("A", "G")``, ``"--"`` -> ``("-", "-")``,
    ``"A"`` -> ``("A", "A")`` (haploid MT/Y).
    """
    if genotype == "--":
        return NO_CALL_MARKER, NO_CALL_MARKER
    if len(genotype) == 2:
        return genotype[0], genotype[1]
    if len(genotype) == 1:
        return genotype, genotype
    logger.warning("Unexpected genotype format %r — treating as no-call", genotype)
    return NO_CALL_MARKER, NO_CALL_MARKER
