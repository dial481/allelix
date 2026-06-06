# Allelix

Open-source command-line toolkit for analyzing raw genotype files from consumer DNA testing services. Format-agnostic ingestion, database-agnostic annotation, offline-first.

> **Status:** Pre-release — six parser formats, four annotators (ClinVar +
> PharmGKB + GWAS Catalog + SNPedia), dual-build ClinVar caches
> (GRCh37 + GRCh38), HTML/JSON/terminal reports, methylation +
> pharmacogenomics focused commands, report diffing. Build
> auto-detection from position data (ADR-0021). No regex on prose
> anywhere in production. Release notes: [`CHANGELOG.md`](CHANGELOG.md).

## Quickstart

Requires Python 3.11+.

```bash
git clone https://github.com/dial481/allelix
cd allelix
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"

# Generate a synthetic test fixture (real fixtures contain personal data; we use mock data)
python tests/generate_mock_data.py

# Show summary statistics for a genotype file
allelix stats tests/fixtures/mock_myhappygenes.txt

# Download reference databases (ClinVar VCF ~100MB compressed; PharmGKB clinical
# annotations ~1MB). Re-runs skip downloads when local + remote signals match;
# use --force to refresh unconditionally.
allelix db update
allelix db status   # see what's cached

# Analyze a genotype file against all ready databases
allelix analyze tests/fixtures/mock_myhappygenes.txt --min-magnitude 5

# Same data, focused subsets
allelix methylation tests/fixtures/mock_myhappygenes.txt
allelix pharmacogenomics tests/fixtures/mock_myhappygenes.txt

# Output to a self-contained HTML or JSON report
allelix analyze tests/fixtures/mock_myhappygenes.txt --output report.html
allelix analyze tests/fixtures/mock_myhappygenes.txt --output report.json
```

## Supported Formats

| Format | Status | Notes |
|---|---|---|
| MyHappyGenes (Tempus) | ✓ | Tab-delimited, 5 columns. **Build is auto-detected** — real-world MHG exports mislabel the header as "build 37.1" while shipping GRCh38 coordinates. Allelix detects from position data and warns on header/data disagreement (ADR-0021). |
| 23andMe | ✓ | Tab-delimited, 4 columns, concatenated genotype. Supports build 36/37/38 from header. I-prefixed probe IDs passed through. |
| AncestryDNA | ✓ | Tab-delimited, 5 columns. Chromosome mapping: 23→X, 24→Y, 25→X (PAR), 26→MT. V1 and V2 chip layouts. |
| Family Tree DNA | ✓ | CSV, double-quoted fields, concatenated genotype. Build 37 default. |
| MyHeritage DNA | ✓ | CSV, same structure as FTDNA. Detected by "MyHeritage" in comment header. Handles double-double-quoted field variant. |
| Living DNA | ✓ | Tab-delimited despite `.csv` extension. Handles AX-, AFFX-prefixed and CHR:POS positional SNP IDs. |

Adding a new format means adding one file to `allelix/parsers/` and registering an instance in the `PARSERS` list in `allelix/parsers/__init__.py`.

### v2 roadmap

| Format | Notes |
|---|---|
| VCF | REF/ALT encoding, `0/1` genotype notation, absence-means-reference semantics. Architecturally different from array parsers — 4-6M variants per file, streaming + batch SQL required. |
| Genome Watchtower | Real-time variant monitoring via database delta feeds. Privacy-preserving: server publishes universal feed, matching happens locally against your deviation set. Replaces full re-analysis with millisecond set intersection. |

## Supported Databases

| Database | Status | Notes |
|---|---|---|
| ClinVar (GRCh37 + GRCh38) | ✓ | Public domain (NCBI). SNVs + indels + multi-allelic sites. **Both builds cached**; `analyze` dispatches by detected build (ADR-0021). Carrier rule (ADR-0007) requires the user to carry the ALT allele. Indel-anchor protection (ADR-0011) prevents single-base array readouts from matching anchor-base indels. |
| PharmGKB | ✓ | CC BY-SA 4.0. Clinical annotations only — single-rsid SNVs; star alleles and haplotypes deferred (ADR-0009). **Primary non-finding filter is the ClinVar REF carrier rule (ADR-0023):** if ClinVar publishes a single-base REF for the rsid and the user is homozygous for it, the row is suppressed. CPIC's `(rsid, base) → function_class` join (ADR-0020) survives as a secondary tier for rsids ClinVar doesn't catalog. Earlier prose tiers (ADR-0013, ADR-0017, ADR-0018) are superseded. |
| CPIC (per-allele function table) | ✓ | Internal data source for the PharmGKB filter. Fetched from `api.cpicpgx.org` at `db update` time. Used to populate the `pharmgkb_allele_function` table — not surfaced to end users as its own annotator. |
| SNPedia | ✓ | CC BY-NC-SA 3.0 US. **Optional — requires a one-time download** via `python scripts/scrape_snpedia.py`. Scrapes both `Category:Is_a_snp` (111,726 pages) and `Category:Is_a_genotype` (104,806 pages) from the MediaWiki API. Stores raw wiki markup; the annotator parses structured genotype templates at query time. If the SNPedia database is absent, analysis runs without it. For commercial use, pass `--exclude-snpedia` or skip the scrape step — either way, `analyze` runs using all other databases and omits SNPedia annotations. |
| GWAS Catalog | ✓ | Public domain (EBI/NHGRI). Trait–SNP associations with p-values and effect sizes. Carrier rule (ADR-0007) requires the user to carry the risk allele. P-value magnitude scoring (ADR-0024) maps continuous p-values to the 0–10 scale; unknown-risk-allele entries fire on rsID match alone but are capped at 3.0. |

### Known PharmGKB limitation: reference-genotype rows where ClinVar and CPIC both lack data

ADR-0022 + ADR-0023: a tiny residual of PharmGKB rows may appear in reports even when the user is homozygous reference. PharmGKB publishes one annotation per genotype including the reference homozygote, and for the reference-homozygote row to be suppressed Allelix needs structured data on the variant from either:

- **ClinVar's REF allele** (the primary filter — see ADR-0023). Covers any rsID ClinVar catalogs.
- **CPIC's per-allele function table** (the secondary fallback — see ADR-0020). Covers rsIDs CPIC has classified.

For the rare rsID where PharmGKB has an annotation but *neither* ClinVar nor CPIC has data, the row emits. These are identifiable by a homozygous-reference genotype combined with "decreased risk," "may have a typical response," or similar comparative language. They are an upstream data gap, not an Allelix bug — we surface them honestly rather than hide them behind a curated exclusion list (which would recreate the maintenance trap the v0.5–v0.7 prose filters were trying to escape).

The CFTR × ivacaftor leak (~30+ rows on real data, pre-v0.7.3) is fixed by the ADR-0023 ClinVar REF check: CPIC's CFTR vocabulary (`"ivacaftor responsive"`) doesn't match the four-class enum the secondary tier expects, but ClinVar publishes REF for every CFTR rsID, so the primary tier catches them universally.

### Known ClinVar upstream data quality issues

Two ClinVar rows in real-world reports are known upstream artifacts, not Allelix bugs:

- **PKD1 rs199476100 GG (Pathogenic/Likely pathogenic, magnitude 8.5).** This is a stop-gained variant with a gnomAD frequency of 0.0005% (7 observations in 1.38 million chromosomes). Homozygosity for this variant is biologically implausible — PKD1 is autosomal dominant and the nonsense variant would be embryonic-lethal or devastating in homozygous state. The chip genotyping call is almost certainly a probe artifact. The code correctly reports what ClinVar says and what the chip reads; the error is upstream of Allelix. Future work: population-frequency filtering could flag ultra-rare variants where the chip call is likely unreliable.

- **IL10 rs1800896 CT (Pathogenic, magnitude 9.0).** This is a common polymorphism (MAF ~20–40%) in the IL-10 promoter. ClinVar's Pathogenic classification comes from a single submitter for hepatitis C susceptibility; a second submitter classifies the same allele as "Uncertain risk allele" for leprosy susceptibility. The ClinVar VCF aggregates across conditions, so the report may pair the Pathogenic classification with the wrong condition. Future work: ClinVar review-status weighting (number of submitters, star rating) could down-weight single-submitter classifications on common variants.

Neither issue affects Allelix's filter logic. Both are inherent to ClinVar's aggregation model and the limitations of array-based genotyping chips.

## Regulatory Posture

Allelix is an informational research tool. It reports classifications made by external databases. It does not independently classify variants, diagnose conditions, or make health recommendations. All variant significance is attributed to its source — Allelix says "ClinVar classifies this variant as pathogenic," never "this variant is pathogenic."

This is not a disclaimer afterthought. It is a design constraint that affects model naming, report wording, and category labeling throughout the codebase.

## Privacy

- No data leaves your machine. No telemetry. No uploads. No analytics.
- Reference databases are downloaded via `allelix db update` and cached locally.
- Analysis runs entirely offline.

## Data Sources & Licensing

Allelix source code is licensed under the **GNU Affero General Public License v3.0 or later** (AGPL-3.0-or-later). Allelix ships with **zero third-party data**. All reference databases are downloaded by the user at runtime via `allelix db update`. Each database retains its original license on the user's machine:

| Database | Source | License | Usage |
|---|---|---|---|
| ClinVar | NCBI | Public domain | No restrictions |
| GWAS Catalog | EBI/NHGRI | Public domain | No restrictions |
| PharmGKB | pharmgkb.org | CC BY-SA 4.0 | Attribution required |
| CPIC | cpicpgx.org | CC BY-SA 4.0 | Attribution required. Per-allele function data fetched from `api.cpicpgx.org` at `db update` time; used internally for the PharmGKB non-finding filter (ADR-0020), not surfaced as its own annotator. |
| SNPedia | snpedia.com | CC BY-NC-SA 3.0 US | Attribution required, **non-commercial only**. Use `--exclude-snpedia` to omit. |

**Commercial users:** SNPedia content is non-commercial. Pass `--exclude-snpedia` to any analysis command, or skip the `scripts/scrape_snpedia.py` step entirely — either way, `analyze` runs using all other databases and omits SNPedia annotations automatically.

### SNPedia data download

SNPedia data is not downloaded by `allelix db update` — it requires a separate one-time scrape:

```bash
python scripts/scrape_snpedia.py
```

This downloads all 216,532 pages (111,726 SNP pages + 104,806 genotype pages) from bots.snpedia.com into `~/.local/share/allelix/snpedia.sqlite` (or `$ALLELIX_DATA_DIR`). The scrape takes 1–4 hours depending on server load. It is resumable — if interrupted, run again to continue. SNPedia is frozen (no new edits since mid-2023), so this is a one-time operation.

If the SNPedia database is not present, `allelix analyze` runs normally using all other databases and prints a note that SNPedia data is not available.

Credit: [jaykobdetar/SNPedia-Scraper](https://github.com/jaykobdetar/SNPedia-Scraper) demonstrated the correct MediaWiki `categorymembers` API approach and published a [Zenodo archive](https://zenodo.org/records/16053572) of the SNP pages. Our scraper extends this by also downloading the 104,806 genotype pages (`Category:Is_a_genotype`), which contain the per-genotype magnitude, repute, and summary data needed for annotation.

### Known SNPedia source data quality notes

SNPedia appears frozen — no edits have been observed since mid-2023. The data below reflects the state of the wiki at scrape time (May 2026) and is unlikely to change.

Of the 104,806 genotype pages in the archive:

- **103 pages have empty or missing allele fields.** These are incomplete entries on the source wiki — the `{{Genotype}}` template was created but the `allele1`/`allele2` fields were never filled in (e.g., `Rs1131692198(;)` with `|allele1=\n|allele2=\n`). All 103 were verified against the live site on 2026-05-21; every one matches the source exactly. The annotator silently skips these — they cannot match any user genotype.

- **1 page has no `{{Genotype}}` template at all.** `Rs1799853(T)` is a malformed single-allele page (`{{is a|genotype}}` instead of a proper genotype template). Skipped by the parser.

- **2 pages have a space before the parenthesis in the title** (`Rs52820871 (G;G)` and `Rs52820871 (G;T)` instead of the standard `Rs52820871(G;G)` format). The annotator handles both title styles.

None of these are scraping errors. They are editorial inconsistencies on the source wiki. The annotator handles all of them correctly: incomplete entries are skipped, variant title formats are matched, and no false annotations are produced.

## Architecture & Design Decisions

The "why" behind major design choices lives in [`docs/adr/`](docs/adr/README.md) as Architecture Decision Records. Read these before proposing changes that touch the parser/annotator interfaces, the regulatory posture, or the data-handling model.

Notable load-bearing ADRs:

- **ADR-0016 — Data Classification Principle.** Classification reads structured fields only. Regex on prose is forbidden in production code.
- **ADR-0020 — CPIC API as the per-allele function source.** The PharmGKB non-finding filter is a table join keyed on `(rsid, base) → clinicalfunctionalstatus`, sourced from CPIC's structured API. Supersedes the prose-extraction tiers from earlier versions (ADR-0017, ADR-0018).
- **ADR-0007 — Genotype matching requires the user to carry the ALT allele.** Applies to ClinVar.
- **ADR-0009 — PharmGKB matches the user's exact normalized diploid call.**
- **ADR-0015 — Mock data generators are the contract.** Fixture shape must mirror real data shape; invariants tested.

Release history: see [`CHANGELOG.md`](CHANGELOG.md).

## Development

```bash
source .venv/bin/activate
pip install -e ".[dev]"

# One-time: install git hooks that block bad commits and pushes
pre-commit install --hook-type pre-commit --hook-type pre-push

ruff check .
ruff format --check .
pytest
```

The hooks (`.pre-commit-config.yaml`) enforce:

- **pre-commit**: `ruff check` + `ruff format --check` (fast lint/format gate)
- **pre-push**: `pytest` (full suite + 92% coverage floor)

If a commit or push is blocked, fix the underlying problem rather than skipping the hook. `git commit --no-verify` and `git push --no-verify` exist but should be used only when you know exactly why.

## License

AGPL-3.0-or-later. See `LICENSE`.
