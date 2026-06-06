# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Abstract base class for reference-database annotators."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, ClassVar

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path
    from types import TracebackType

    from allelix.models import Annotation, Variant


def is_clinvar_homref(
    variant: Variant,
    clinvar_ref_provider: Callable[[str, str], str | None] | None,
) -> bool:
    """Return True if the variant is homozygous reference per ClinVar (ADR-0023)."""
    if clinvar_ref_provider is None:
        return False
    ref = clinvar_ref_provider(variant.rsid, variant.build)
    return ref is not None and len(ref) == 1 and variant.allele1 == ref and variant.allele2 == ref


class Annotator(ABC):
    """Base class for all reference database annotators.

    Annotators bind to a `data_dir` at construction. `setup()` performs the
    one-time download/parse into the cache. `is_ready()` reports whether the
    cache exists and is queryable. `annotate(variant)` returns zero or more
    `Annotation` objects for the variant — checking both rsid AND genotype, per
    the regulatory posture (ADR-0003) and the genotype-matching rule (ADR-0007).

    Annotators hold resources (SQLite connections, file handles). Always close
    them via `close()` or the context manager protocol; the CLI uses
    `contextlib.ExitStack` to guarantee deterministic cleanup.

    Attributes:
        name: Lowercase identifier (e.g., "clinvar").
        display_name: Human-readable name ("ClinVar").
        attribution: Display label used in user-facing reports ("ClinVar").
            Equal to `display_name` for first-party single-source annotators.
        requires_download: Whether `setup()` needs network/disk space.
    """

    name: ClassVar[str]
    display_name: ClassVar[str]
    attribution: ClassVar[str]
    requires_download: ClassVar[bool] = True

    def __init__(self, data_dir: Path) -> None:
        """Bind the annotator to a data directory (created elsewhere)."""
        self.data_dir = data_dir

    def __enter__(self) -> Annotator:
        """Return self for `with` usage."""
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> None:
        """Release any open resources via `close()`."""
        self.close()

    @abstractmethod
    def setup(self) -> None:
        """Download and prepare the reference database. Idempotent."""
        ...

    @abstractmethod
    def annotate(self, variant: Variant) -> list[Annotation]:
        """Return all annotations for this variant.

        Implementations MUST verify both rsid AND genotype — presence in the
        database is not enough. The user must carry the flagged allele.
        """
        ...

    @abstractmethod
    def is_ready(self) -> bool:
        """Whether the local cache exists and is queryable."""
        ...

    @abstractmethod
    def version(self) -> str | None:
        """Return the cached database version, or None if not set up."""
        ...

    @abstractmethod
    def close(self) -> None:
        """Release any open resources (database connections, file handles)."""
        ...

    @abstractmethod
    def fetch_remote_signal(self) -> str | None:
        """Fetch a small remote freshness signal (md5 hash, ETag, Last-Modified).

        Implementations MUST return a prefixed, opaque string (e.g.
        `"md5:abcdef…"`, `"lm:Wed, 21 Oct 2025 …"`, `"etag:…"`) so that a
        future server-side switch in signal type triggers a refresh
        rather than a silent miss.

        Returns None on any failure (network error, timeout, missing
        header, source doesn't expose a signal). Never raises — `db update`
        treats None as "can't verify freshness" and falls through to skip
        with a notice. See ADR-0012.
        """
        ...

    @abstractmethod
    def cached_remote_signal(self) -> str | None:
        """Return the remote signal stored at last successful download, or None.

        Returns None if the cache is missing entirely OR if the cache was
        written by a pre-v0.4.2 release that didn't capture signals.
        """
        ...
