# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Interpreter version stamps for annotator cache invalidation.

Increment the constant for an annotator when its emit/suppression logic
changes in a way that should invalidate prior reports built against
existing caches. The stamp is appended to ``database_versions.remote_signal``
as ``|iv:N`` so ``is_ready()`` can reject stale caches without forcing a
full re-download.
"""

CLINVAR_INTERPRETER_VERSION = 1
PHARMGKB_INTERPRETER_VERSION = 1
