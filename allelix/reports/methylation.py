# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Methylation pathway gene panel.

The set is intentionally small and curated — covering the one-carbon /
folate / methylation cycle genes most often discussed in consumer methylation
reports. Not medical guidance; see ADR-0003 (regulatory posture).
"""

from __future__ import annotations

# Folate / one-carbon / methylation cycle genes. Add via PR + ADR if expanding.
METHYLATION_PANEL_GENES: frozenset[str] = frozenset(
    {
        "MTHFR",
        "MTR",
        "MTRR",
        "COMT",
        "CBS",
        "BHMT",
        "BHMT2",
        "MTHFD1",
        "MTHFD1L",
        "AHCY",
        "MAT1A",
        "MAT2A",
        "MAT2B",
        "TCN1",
        "TCN2",
        "FUT2",
        "GSTM1",
        "GSTP1",
        "DNMT1",
        "DNMT3A",
        "DNMT3B",
        "PEMT",
        "FOLR1",
        "FOLR2",
        "SLC19A1",
        "SHMT1",
        "SHMT2",
    }
)
