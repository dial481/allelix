# Full Test Protocol

External reviewer checklist for verifying an allelix release against real data.

**Requirements:** Fast machine, fast internet, ~50 GB free disk space.
Estimated wall-clock time: 30–45 minutes (most of it is database downloads).

## 1. Environment setup

```bash
git clone git@github.com:dial481/allelix.git
cd allelix
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Verify Python 3.11+:

```bash
python --version
```

## 2. Unit test suite (synthetic data)

Run the full test suite against synthetic fixtures. No network access
required — all mock data is committed.

```bash
python -m pytest tests/ -x --tb=short
```

**Expected:** 1090 tests pass, 0 failures.

Check lint:

```bash
ruff check . && ruff format --check .
```

**Expected:** All checks passed, 0 files reformatted.

## 3. Download all databases

```bash
allelix db update
```

This downloads ClinVar (GRCh37 + GRCh38), PharmGKB, GWAS Catalog,
gnomAD (~6 GB), AlphaMissense (~8 GB), and SNPedia from HuggingFace.

**Expected:** All 6 annotators show green checkmarks. No errors.

Verify status:

```bash
allelix db status
```

**Expected:** All annotators show "yes" in the Ready column with
version strings and record counts. SNPedia should show ~104K records.

## 4. Fetch real test data

```bash
bash scripts/fetch_testdata.sh
```

This downloads CC0 public-domain genotype files from the GitHub release
and the GWAS Catalog from EBI (~66 MB).

**Expected:** `test_data/real/` and `test_data/transcoded/` directories
populated. GWAS catalog zip present at `test_data/gwas_catalog.zip`.

## 5. Analyze real genotype files

Run analysis on each format against the live-downloaded databases.

### 5a. 23andMe

```bash
allelix analyze test_data/real/23andme/user1190_v5.txt --output /tmp/allelix-review/user1190_23andme.json
```

**Expected:** Exit code 0. JSON report written. Should contain ClinVar,
PharmGKB, GWAS, SNPedia, gnomAD, and AlphaMissense annotations. Check
that annotation count is in the hundreds (varies by database version).

### 5b. MHG / Tempus

```bash
allelix analyze test_data/real/mhg/user1190.txt --output /tmp/allelix-review/user1190_mhg.json
```

**Expected:** Exit code 0. JSON report written. This file is a clean
GRCh37 transcode of user1190_v5.txt — no build mismatch expected.
(The mismatch fixture is `edge_cases/mhg_grch38_with_grch37_header.txt`,
tested in step 14.)

### 5c. AncestryDNA

```bash
# Pick any one file from the directory
allelix analyze "$(find test_data/real/ancestrydna -maxdepth 1 -type f | head -1)" \
  --output /tmp/allelix-review/ancestrydna.json
```

**Expected:** Exit code 0. JSON report written.

### 5d. FTDNA

```bash
# Pick any one file from the directory
allelix analyze "$(find test_data/real/ftdna -maxdepth 1 -type f | head -1)" \
  --output /tmp/allelix-review/ftdna.json
```

### 5e. Living DNA

```bash
allelix analyze test_data/real/livingdna/user1190.csv --output /tmp/allelix-review/user1190_livingdna.json
```

### 5f. MyHeritage

```bash
allelix analyze test_data/real/myheritage/user1190.csv --output /tmp/allelix-review/user1190_myheritage.json
```

## 6. Cross-parser identity check

The user1190 genotype exists in 6 format representations. All should
produce identical annotation sets (same rsIDs, same significance, same
sources). The exact annotation count depends on database versions, but
the counts must match across formats.

```bash
mkdir -p /tmp/allelix-review
for f in \
  test_data/real/23andme/user1190_v5.txt \
  test_data/real/mhg/user1190.txt \
  test_data/real/livingdna/user1190.csv \
  test_data/real/myheritage/user1190.csv \
  test_data/transcoded/user1190_as_ancestrydna.txt \
  test_data/transcoded/user1190_as_ftdna.csv; do
  echo "=== $f ==="
  allelix analyze "$f" --exclude-snpedia --output /tmp/allelix-review/$(basename "$f").json 2>&1 | tail -3
done
```

Then compare annotation counts:

```bash
for f in /tmp/allelix-review/user1190_*.json; do
  echo "$(basename $f): $(python3 -c "import json; print(len(json.load(open('$f'))['annotations']))")"
done
```

**Expected:** All 6 files produce the same annotation count. Any
discrepancy is a parser or build-detection bug.

## 7. Multi-allelic enrichment accuracy (issue #25)

Verify that enrichment lookups use exact alt-allele matching, not
MAX-aggregated fallback.

```bash
allelix analyze test_data/real/23andme/user1190_v5.txt --output /tmp/allelix-review/enrichment_check.json
python3 -c "
import json
data = json.load(open('/tmp/allelix-review/enrichment_check.json'))
for a in data['annotations']:
    if a.get('am_pathogenicity') is not None and a.get('alt'):
        print(f\"{a['rsid']} alt={a['alt']} am={a['am_pathogenicity']:.3f} {a['am_class']}\")
" | head -20
```

**Expected:** AM scores correspond to the user's specific alt allele,
not the site-wide MAX. Spot-check a few rsIDs against the AlphaMissense
source data if available.

## 8. Report formats

### 8a. HTML report

```bash
allelix analyze test_data/real/23andme/user1190_v5.txt --output /tmp/allelix-review/report.html
```

Open `/tmp/allelix-review/report.html` in a browser. Verify:

- Table renders without horizontal overflow
- rsID column is sticky when scrolling
- Columns are sortable (click headers)
- Review Status column appears for ClinVar rows
- Pop. Freq column shows gnomAD frequencies
- AM Score column shows AlphaMissense scores
- PharmGKB AM scores show dimmed caveat indicator
- Row borders are color-coded (red = pathogenic, green = benign)
- "Reading This Report" section is present
- Regulatory notice is present

### 8b. Terminal report

```bash
allelix analyze test_data/real/23andme/user1190_v5.txt 2>&1 | head -50
```

**Expected:** Rich-formatted table with colored output. All columns
present.

### 8c. JSON report

```bash
python3 -c "
import json, sys
data = json.load(open('/tmp/allelix-review/enrichment_check.json'))
print(f\"Schema version: {data.get('schema_version')}\")
print(f\"Annotations: {len(data['annotations'])}\")
print(f\"Sources: {set(a['source'] for a in data['annotations'])}\")
has_af = sum(1 for a in data['annotations'] if a.get('allele_frequency') is not None)
has_am = sum(1 for a in data['annotations'] if a.get('am_pathogenicity') is not None)
print(f\"With gnomAD freq: {has_af}\")
print(f\"With AM score: {has_am}\")
"
```

**Expected:** Schema version 3. Multiple sources present. gnomAD and
AM enrichment counts > 0.

## 9. Stats, extract, and focused reports

```bash
allelix stats test_data/real/23andme/user1190_v5.txt
allelix extract --snps rs1801133,rs429358,rs7412 test_data/real/23andme/user1190_v5.txt
```

**Expected:** Stats shows SNP count, no-call rate, het rate. Extract
returns the requested SNPs with genotypes.

### 9a. Focused subcommands

```bash
allelix methylation test_data/real/23andme/user1190_v5.txt
```

**Expected:** Methylation pathway report with annotations from the
methylation gene panel. Non-zero annotation count.

```bash
allelix pharmacogenomics test_data/real/23andme/user1190_v5.txt
```

**Expected:** PharmGKB-focused report. Non-zero annotation count.

### 9b. Compare

```bash
allelix compare test_data/real/23andme/user1190_v5.txt test_data/real/myheritage/user1190.csv
```

**Expected:** Per-chromosome concordance table. Coverage overlap stats.
High concordance expected (same biology, different format).

## 10. Config system

```bash
allelix config show
allelix config set license.commercial true
allelix config show
allelix analyze test_data/real/23andme/user1190_v5.txt 2>&1 | grep -i "snpedia\|skipping"
allelix config set license.commercial false
allelix config show
```

**Expected:** With `license.commercial = true`, SNPedia is excluded
from analysis automatically. After setting back to `false`, SNPedia
is included again.

## 11. Diff command

```bash
allelix analyze test_data/real/23andme/user1190_v5.txt --output /tmp/allelix-review/baseline.json
allelix analyze test_data/real/23andme/user1190_v5.txt --output /tmp/allelix-review/current.json --diff /tmp/allelix-review/baseline.json
```

**Expected:** Diff reports no changes (same input, same databases).

## 12. Database update signals

```bash
allelix db update
```

**Expected:** Most annotators show "already current". Per-annotator
states:

- ClinVar, GWAS Catalog (server-driven): "already current" or "can't
  be verified" (ETag/sidecar-dependent)
- PharmGKB (server-driven, CPIC-API dependent): "already current" or
  "can't be verified"
- gnomAD, AlphaMissense, SNPedia (code-driven, ADR-0030): always
  "already current" — refresh only via `--force` or code bump of
  pinned commit SHA

No re-downloads.

```bash
allelix db update --force
```

**Expected:** All annotators re-download and show green checkmarks.
Note: `--force` semantics differ by tier. Server-driven sources
override a "signal matches" skip; code-driven sources have no
signal-match path to override — `--force` is the only way to
re-trigger their download because pinned URLs are deterministic.
See ADR-0030.

## 13. GWAS Catalog real-data sanity (slow tests)

These tests load the real 795K-record GWAS Catalog and verify that the
magnitude scoring formula produces bounded output.

```bash
python -m pytest tests/test_end_to_end.py -k "TestRealDataGwasSanity" -v
```

**Expected:** 2 tests pass. Default floor (9.0) keeps output under 50
rows. Old floor (7.0) produces more output than new floor.

## 14. Edge case files

```bash
# Build mismatch detection (analyze runs the build-detection pipeline; stats does not)
allelix analyze test_data/edge_cases/mhg_grch38_with_grch37_header.txt 2>&1 | grep -i "mismatch\|build"
# Expected: Build mismatch warning (header claims GRCh37, positions are GRCh38)

# P-A: canonical header tightening — this file should NOT be recognized as 23andMe
allelix stats test_data/edge_cases/23andme_lookalike_rejected_by_PA.txt 2>&1
# Expected: "No parser recognized" error

# Genes for Good — 23andMe-format export from a different service
allelix stats test_data/edge_cases/23andme_format_from_genes_for_good_service.txt
# Expected: Recognized as 23andMe format, stats displayed

# GRCh36 FTDNA file (analyze detects build from positions; stats shows parser default)
allelix analyze test_data/edge_cases/ftdna_grch36_positions.csv 2>&1 | grep -i "grch36\|build"
# Expected: GRCh36 detected. ClinVar skipped (no GRCh36 cache).

# Unsupported formats
allelix stats test_data/edge_cases/unsupported_decodeme.txt 2>&1
allelix stats test_data/edge_cases/unsupported_23andme_exome_vcf.txt 2>&1
# Expected: "No parser recognized" for both
```

## 15. Cleanup

```bash
rm -rf /tmp/allelix-review
```

Optionally remove downloaded databases to free ~15 GB:

```bash
rm -rf ~/.local/share/allelix/
```

## Pass criteria

All of the following must be true:

- [ ] Unit test suite: 1090 passed, 0 failed
- [ ] Ruff lint + format: zero warnings
- [ ] `db update` downloads all 6 annotators without errors
- [ ] `db status` shows all annotators ready with version and record count
- [ ] All 6 parser formats produce successful analysis
- [ ] Cross-parser identity: same annotation count across all user1190 representations
- [ ] HTML report renders correctly in a browser
- [ ] JSON report has schema version 3 with gnomAD + AM enrichment
- [ ] Config system correctly gates SNPedia on `license.commercial`
- [ ] Edge case files produce expected behavior
- [ ] `db update` (second run) skips already-current databases
- [ ] GWAS Catalog slow tests pass
- [ ] `methylation`, `pharmacogenomics`, `compare` subcommands produce output
