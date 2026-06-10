# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Tests for the ClinVar annotator."""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename, clinvar_record_name
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.manager import load_clinvar_vcf
from allelix.models import Variant

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture
def annotator(clinvar_data_dir: Path):
    """Yield an annotator and ensure its SQLite connection is closed.

    N-1: without explicit teardown, _connection() opens a sqlite3.Connection
    that's only reaped by GC. Yield + close() pins the contract every test
    relies on.
    """
    ann = ClinVarAnnotator(clinvar_data_dir)
    try:
        yield ann
    finally:
        ann.close()


@pytest.fixture
def annotator_with_benign(clinvar_data_dir: Path):
    """Annotator that includes benign/likely_benign annotations."""
    ann = ClinVarAnnotator(clinvar_data_dir, include_benign=True)
    try:
        yield ann
    finally:
        ann.close()


class TestSetupAndStatus:
    def test_unconfigured_is_not_ready(self, tmp_path: Path):
        ann = ClinVarAnnotator(tmp_path)
        assert ann.is_ready() is False
        assert ann.version() is None

    def test_configured_is_ready(self, annotator: ClinVarAnnotator):
        assert annotator.is_ready() is True
        assert annotator.version() is not None


class TestSignalGuard:
    def test_setup_aborts_when_signal_fetch_fails(self, tmp_path: Path, monkeypatch):
        """setup() raises RuntimeError when remote signal is None."""
        ann = ClinVarAnnotator(tmp_path)
        monkeypatch.setattr(
            ClinVarAnnotator, "_fetch_remote_signal_for", staticmethod(lambda _build: None)
        )
        with pytest.raises(RuntimeError, match="cannot verify remote freshness signal"):
            ann.setup()


class TestInterpreterVersionStamp:
    """CLINVAR_INTERPRETER_VERSION stamp in cache's local_version_tag."""

    def test_is_ready_accepts_matching_iv_stamp(self, annotator: ClinVarAnnotator):
        """Freshly loaded cache has the current iv stamp — is_ready returns True."""
        assert annotator.is_ready() is True

    def test_is_ready_rejects_cache_without_tag(
        self, tmp_path: Path, mock_clinvar_grch37_vcf: Path
    ):
        """Cache with no local_version_tag is self-healed by one-shot migration."""
        build = "GRCh37"
        db_path = tmp_path / clinvar_db_filename(build)
        load_clinvar_vcf(
            mock_clinvar_grch37_vcf,
            db_path,
            source_url="test://mock",
            record_name=clinvar_record_name(build),
        )
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET remote_signal = 'md5:abc', "
                "local_version_tag = NULL WHERE name = ?",
                (clinvar_record_name(build),),
            )
            conn.commit()
        ann = ClinVarAnnotator(tmp_path, builds=(build,))
        assert ann.is_ready() is True
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            tag = conn.execute(
                "SELECT local_version_tag FROM database_versions WHERE name = ?",
                (clinvar_record_name(build),),
            ).fetchone()[0]
            assert tag == f"iv:{CLINVAR_INTERPRETER_VERSION}"

    def test_is_ready_rejects_old_iv_stamp(self, tmp_path: Path, mock_clinvar_grch37_vcf: Path):
        """Cache stamped with an older iv version is rejected."""
        build = "GRCh37"
        db_path = tmp_path / clinvar_db_filename(build)
        load_clinvar_vcf(
            mock_clinvar_grch37_vcf,
            db_path,
            source_url="test://mock",
            record_name=clinvar_record_name(build),
        )
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET local_version_tag = 'iv:0' WHERE name = ?",
                (clinvar_record_name(build),),
            )
            conn.commit()
        ann = ClinVarAnnotator(tmp_path, builds=(build,))
        assert ann.is_ready() is False


class TestGenotypeMatching:
    """ADR-0007: ClinVar entries trigger only when the user carries ALT."""

    def test_heterozygous_carrier_triggers(self, annotator: ClinVarAnnotator):
        # mock ClinVar: rs1801133, REF=G, ALT=A, Pathogenic
        v = Variant("rs1801133", "1", 11796321, "G", "A")
        results = annotator.annotate(v)
        assert len(results) == 1
        a = results[0]
        assert a.significance == "clinvar_pathogenic"
        assert a.attribution == "ClinVar"
        assert a.source == "clinvar"
        assert a.category == "clinical"
        assert a.gene == "MTHFR"
        # ADR-0023: genotype_match shows the user's diploid (sorted), not
        # the matched ALT base. G/A → "AG".
        assert a.genotype_match == "AG"
        assert a.magnitude == 9.0

    def test_homozygous_alt_triggers(self, annotator: ClinVarAnnotator):
        # mock ClinVar: rs4680, REF=G, ALT=A, Drug_response
        v = Variant("rs4680", "22", 19963748, "A", "A")
        results = annotator.annotate(v)
        assert len(results) == 1
        assert results[0].significance == "clinvar_drug_response"
        assert results[0].magnitude == 6.5

    def test_homozygous_reference_does_not_trigger(self, annotator: ClinVarAnnotator):
        # mock ClinVar: rs121918506, REF=G, ALT=T, Pathogenic; mock has G/G
        v = Variant("rs121918506", "17", 7577538, "G", "G")
        assert annotator.annotate(v) == []

    def test_no_call_does_not_trigger(self, annotator: ClinVarAnnotator):
        v = Variant("rs1801133", "1", 11796321, "-", "-")
        assert annotator.annotate(v) == []

    def test_asymmetric_no_call_does_not_trigger(self, annotator: ClinVarAnnotator):
        """r-2: one good allele + one no-call must short-circuit before lookup.

        Catches mutations like `if variant.is_no_call` → `if variant.allele1 == "-"`
        that pass the both-no-call test but leak through here.
        """
        v_left = Variant("rs1801133", "1", 11796321, "-", "A")
        v_right = Variant("rs1801133", "1", 11796321, "A", "-")
        assert annotator.annotate(v_left) == []
        assert annotator.annotate(v_right) == []

    def test_unknown_rsid_does_not_trigger(self, annotator: ClinVarAnnotator):
        v = Variant("rs999000111", "1", 1000, "A", "T")
        assert annotator.annotate(v) == []


class TestAttribution:
    """ADR-0003: Significance and attribution must be source-prefixed."""

    def test_all_annotations_attribute_to_clinvar(self, annotator: ClinVarAnnotator):
        v = Variant("rs1801133", "1", 11796321, "G", "A")
        results = annotator.annotate(v)
        for a in results:
            assert a.attribution == "ClinVar"
            assert a.significance.startswith("clinvar_")
            assert a.category == "clinical"

    def test_description_attributes_to_clinvar(self, annotator: ClinVarAnnotator):
        v = Variant("rs1801133", "1", 11796321, "G", "A")
        results = annotator.annotate(v)
        assert results[0].description.startswith("ClinVar classifies")

    def test_review_status_populated(self, annotator: ClinVarAnnotator):
        """CLNREVSTAT is surfaced on the Annotation."""
        v = Variant("rs1801133", "1", 11796321, "G", "A")
        results = annotator.annotate(v)
        assert len(results) >= 1
        assert results[0].review_status == "criteria_provided,_single_submitter"


class TestRegistryMetadata:
    def test_class_attributes(self):
        assert ClinVarAnnotator.name == "clinvar"
        assert ClinVarAnnotator.display_name == "ClinVar"
        assert ClinVarAnnotator.attribution == "ClinVar"
        assert ClinVarAnnotator.requires_download is True


class TestIndelMatching:
    """M-4: ClinVar contains pathogenic indels (e.g., CFTR ΔF508). Must match."""

    def test_indel_carrier_triggers(self, annotator: ClinVarAnnotator):
        # mock ClinVar: rs113993960 REF=CTT ALT=C, Pathogenic CFTR
        v = Variant("rs113993960", "7", 117199644, "CTT", "C")
        results = annotator.annotate(v)
        assert len(results) == 1
        assert results[0].gene == "CFTR"
        # ADR-0023: indel diploid passes through as `"CTT/C"` to keep
        # multi-base alleles readable rather than concatenating them.
        assert results[0].genotype_match == "CTT/C"
        assert results[0].significance == "clinvar_pathogenic"

    def test_indel_homozygous_reference_does_not_trigger(self, annotator: ClinVarAnnotator):
        v = Variant("rs113993960", "7", 117199644, "CTT", "CTT")
        assert annotator.annotate(v) == []


class TestIndelAnchorProtection:
    """ADR-0011: indel rows must NOT fire on single-base array readouts.

    ClinVar encodes indels with anchor-base notation (REF=CTT ALT=C). Array
    parsers report single bases at probe positions. Pre-v0.4.2 the carrier
    rule's `alt in {allele1, allele2}` matched ClinVar's single-character
    anchor against an array's single-character readout, producing categorical
    false-positive "Pathogenic" calls in cancer-predisposition genes for
    users who carried only the wild-type sequence.
    """

    def test_array_single_base_does_not_fire_on_indel_row(self, annotator: ClinVarAnnotator):
        # mock fixture: rs113993960 REF=CTT ALT=C Pathogenic CFTR.
        # Array reads a single C at the probe position; the user does NOT
        # carry the deletion. Pre-v0.4.2 incorrectly fired.
        v = Variant("rs113993960", "7", 117199644, "C", "C")
        assert annotator.annotate(v) == []

    def test_indel_calling_parser_still_fires(self, annotator: ClinVarAnnotator):
        # A multi-base genotype like CTT/C indicates a parser that actually
        # calls indels (future VCF parser). Indel matching must still work.
        v = Variant("rs113993960", "7", 117199644, "CTT", "C")
        results = annotator.annotate(v)
        assert len(results) == 1
        assert results[0].significance == "clinvar_pathogenic"

    def test_homozygous_alt_indel_still_fires_for_multibase_parser(
        self, annotator: ClinVarAnnotator
    ):
        # Hypothetical homozygous deletion. The user's genotype carries the
        # multi-base form on at least one side, so the indel filter doesn't
        # short-circuit; the carrier rule still applies.
        v = Variant("rs113993960", "7", 117199644, "C", "CTT")
        results = annotator.annotate(v)
        assert len(results) == 1


class TestMultiAllelicMatching:
    """C-2: Multi-allelic ClinVar rows must match per-ALT, not as the joined string."""

    def test_carrier_of_pathogenic_alt_triggers(self, annotator: ClinVarAnnotator):
        # mock ClinVar: rs1065852 ALT=A,C with CLNSIG=Drug_response|Benign.
        # MHG fixture has G/A — carries A only.
        v = Variant("rs1065852", "22", 42526694, "G", "A")
        results = annotator.annotate(v)
        # Should match exactly the A-allele record (Drug_response), not the C one.
        sigs = {r.significance for r in results}
        assert "clinvar_drug_response" in sigs
        assert "clinvar_benign" not in sigs

    def test_carrier_of_benign_alt_only(self, annotator_with_benign: ClinVarAnnotator):
        # User carries G/C — only the C-allele record should fire (Benign).
        v = Variant("rs1065852", "22", 42526694, "G", "C")
        results = annotator_with_benign.annotate(v)
        sigs = {r.significance for r in results}
        assert "clinvar_benign" in sigs
        assert "clinvar_drug_response" not in sigs


class TestRemoteSignal:
    """Freshness signal: ClinVar uses the .md5 sidecar file (ADR-0012)."""

    def test_fetch_returns_md5_prefixed_signal(self, annotator: ClinVarAnnotator, monkeypatch):
        """ADR-0021: signal is composite across managed builds."""
        from allelix.annotators import clinvar as clinvar_module

        monkeypatch.setattr(
            clinvar_module,
            "fetch_remote_text",
            lambda url: "abcdef0123456789  clinvar.vcf.gz\n",
        )
        # Default annotator manages both builds → composite signal.
        assert annotator.fetch_remote_signal() == (
            "GRCh37:md5:abcdef0123456789|GRCh38:md5:abcdef0123456789"
        )

    def test_fetch_returns_none_on_network_error(self, annotator: ClinVarAnnotator, monkeypatch):
        from allelix.annotators import clinvar as clinvar_module

        monkeypatch.setattr(clinvar_module, "fetch_remote_text", lambda url: None)
        assert annotator.fetch_remote_signal() is None

    def test_fetch_returns_none_on_empty_md5_body(self, annotator: ClinVarAnnotator, monkeypatch):
        from allelix.annotators import clinvar as clinvar_module

        monkeypatch.setattr(clinvar_module, "fetch_remote_text", lambda url: "   \n")
        assert annotator.fetch_remote_signal() is None

    def test_cached_returns_none_for_unconfigured(self, tmp_path: Path):
        ann = ClinVarAnnotator(tmp_path)
        assert ann.cached_remote_signal() is None

    def test_cached_returns_none_for_v041_cache(self, annotator: ClinVarAnnotator):
        """v0.4.1 caches were populated without a remote_signal column."""
        # The clinvar_data_dir fixture writes via load_clinvar_vcf without
        # passing remote_signal, so the column exists (new schema) but the
        # value is NULL — cached_remote_signal should return None.
        assert annotator.cached_remote_signal() is None

    def test_cached_round_trip_after_setup(self, tmp_path: Path, mock_clinvar_vcf: Path):
        """ADR-0021: composite cached signal is `GRCh37:<sig>|GRCh38:<sig>`.

        For a single-build annotator the composite collapses to one part.
        """
        from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
        from allelix.databases.manager import load_clinvar_vcf

        load_clinvar_vcf(
            mock_clinvar_vcf,
            tmp_path / clinvar_db_filename("GRCh37"),
            source_url="test",
            remote_signal="md5:deadbeef",
            record_name=clinvar_record_name("GRCh37"),
        )
        ann = ClinVarAnnotator(tmp_path, builds=("GRCh37",))
        try:
            assert ann.cached_remote_signal() == "GRCh37:md5:deadbeef"
        finally:
            ann.close()


class TestConstructorValidation:
    def test_unsupported_build_raises(self, tmp_path: Path):
        with pytest.raises(ValueError, match="Unsupported"):
            ClinVarAnnotator(tmp_path, builds=("GRCh99",))


class TestCloseable:
    """C-1: ClinVarAnnotator must release its SQLite connections deterministically."""

    def test_close_releases_connection(self, annotator: ClinVarAnnotator):
        # Touch the connection
        annotator.annotate(Variant("rs1801133", "1", 11796321, "G", "A"))
        assert annotator._conns, "expected at least one open per-build connection"
        annotator.close()
        assert annotator._conns == {}

    def test_close_is_idempotent(self, annotator: ClinVarAnnotator):
        annotator.close()
        annotator.close()  # must not raise

    def test_context_manager_closes_on_exit(self, clinvar_data_dir: Path):
        ann = ClinVarAnnotator(clinvar_data_dir)
        with ann as bound:
            assert bound is ann
            ann.annotate(Variant("rs1801133", "1", 11796321, "G", "A"))
            assert ann._conns
        assert ann._conns == {}
