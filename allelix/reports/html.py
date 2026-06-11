# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 dial481
"""Self-contained HTML report renderer.

The output is a single `.html` file with inline CSS and no external assets,
suitable for emailing or hosting on a static web mount. Per ADR-0003, every
row carries its source attribution and the page header restates the
informational/non-diagnostic posture.
"""

from __future__ import annotations

import html
from datetime import UTC, datetime
from typing import TYPE_CHECKING

from allelix import __version__
from allelix.reports import REGULATORY_NOTICE, atomic_write_text
from allelix.reports._pipeline import rollup_gwas_duplicates

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path

    from allelix.models import Annotation
    from allelix.reports._pipeline import AnalysisResult
    from allelix.reports.diff import DiffResult


_CSS = """
* { box-sizing: border-box; }
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  margin: 2rem auto; max-width: 1200px; padding: 0 1rem; color: #222;
}
h1 { margin-bottom: .25rem; }
.subtitle { color: #666; margin-top: 0; }
.notice {
  background: #fff8e1; border-left: 4px solid #f9a825;
  padding: 1rem; margin: 1.5rem 0; border-radius: 4px; font-size: .95rem;
}
.notice-warn {
  background: #fff3e0; border-left: 4px solid #e65100;
  padding: 1rem; margin: 1.5rem 0; border-radius: 4px; font-size: .95rem;
}
.education {
  background: #f5f5f5; border-left: 4px solid #90a4ae;
  padding: 1rem 1.25rem; margin: 1.5rem 0; border-radius: 4px; font-size: .9rem;
}
.education h2 { margin-top: 0; font-size: 1.05rem; }
.education h3 { font-size: .95rem; margin: .75rem 0 .25rem; }
.education p { margin: .35rem 0; }
details.education summary { cursor: pointer; font-size: 1.05rem; }
details.education[open] summary { margin-bottom: .5rem; }
.summary {
  display: flex; flex-wrap: wrap; gap: .75rem; margin-bottom: 2rem;
}
.card {
  background: #f5f5f5; padding: .75rem 1rem; border-radius: 4px;
  flex: 1 1 150px; min-width: 120px;
}
.card .label {
  font-size: .8rem; color: #666;
  text-transform: uppercase; letter-spacing: .05em;
}
.card .value { font-size: 1.2rem; font-weight: 600; }
.table-wrap { overflow-x: auto; margin-bottom: 2rem; }
table { width: 100%; border-collapse: collapse; font-size: .9rem; }
th, td {
  text-align: left; padding: .55rem .5rem;
  border-bottom: 1px solid #eee; vertical-align: top;
  white-space: nowrap;
}
td.desc-cell {
  white-space: normal; max-width: 400px; word-break: break-word;
}
th {
  background: #fafafa; position: sticky; top: 0; font-weight: 600;
  cursor: pointer; user-select: none; z-index: 2;
}
th:hover { background: #f0f0f0; }
th .sort-arrow { font-size: .7rem; margin-left: .25rem; color: #999; }
td.col-rsid {
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  position: sticky; left: 0; background: #fff; z-index: 1;
}
th:first-child { position: sticky; left: 0; z-index: 3; }
tr:hover td.col-rsid { background: #fcfcfc; }
.gene { font-weight: 600; color: #5e35b1; }
.source {
  font-size: .8rem; padding: .15rem .5rem;
  border-radius: 4px; background: #e3f2fd; color: #1565c0;
}
.bar {
  display: inline-block; height: 6px; background: #1565c0;
  border-radius: 3px; vertical-align: middle; margin-right: .35rem;
}
.condition { color: #555; font-size: .85rem; }
.refs-toggle { font-size: .8rem; color: #888; margin-top: .25rem; }
.refs-toggle summary { cursor: pointer; }
.empty { padding: 2rem; text-align: center; color: #666; font-style: italic; }
tr.repute-bad { border-left: 4px solid #c62828; background: #fce4e4; }
tr.repute-good { border-left: 4px solid #2e7d32; background: #e8f5e9; }
tr.repute-neutral { border-left: 4px solid #bdbdbd; }
tr.repute-bad td.col-rsid { background: #fce4e4; }
tr.repute-good td.col-rsid { background: #e8f5e9; }
.am-pathogenic { color: #c62828; font-weight: 600; }
.am-benign { color: #2e7d32; font-weight: 600; }
.am-ambiguous { color: #f9a825; font-weight: 600; }
.am-score { color: #999; }
.cadd-high { color: #c62828; font-weight: 600; }
.cadd-med { color: #e65100; font-weight: 600; }
.cadd-low { color: #999; }
tr.diff-new { background: #e8f5e9; }
tr.diff-changed { background: #fff3e0; }
tr.diff-removed { background: #fafafa; text-decoration: line-through; color: #999; }
footer {
  margin-top: 3rem; padding-top: 1rem; border-top: 1px solid #eee;
  font-size: .8rem; color: #888;
}
"""

_EDUCATION_SECTION = """
<div class="education">
<h2>Reading This Report</h2>

<h3>Pseudogene cross-hybridization</h3>
<p>Some genes have known pseudogenes with high sequence similarity
(e.g., PKD1 has six pseudogenes sharing &gt;97% identity). Array-based
genotyping probes can cross-hybridize, producing false genotype calls.
If a result seems inconsistent with your health history, confirmatory
testing by a different method is recommended.</p>

<h3>ClinVar aggregation</h3>
<p>ClinVar aggregates submissions from multiple sources. Different submitters
may classify the same variant differently, and significance labels may be
paired with conditions from a different submitter&rsquo;s entry. Single-submitter
entries carry less weight than expert-panel-reviewed classifications.</p>

<h3>Carrier vs. affected</h3>
<p>A variant classified as pathogenic in a recessive condition requires two
copies to cause disease. If you are heterozygous (one copy), you are a carrier.
Carrier status does not typically cause symptoms but may be relevant for family
planning.</p>

<h3>Confirmatory testing</h3>
<p>No genotyping platform is 100% accurate. Clinically significant findings
should be confirmed with an independent method before making medical
decisions.</p>
</div>
"""

_MAGNITUDE_LEGEND = """
<details class="education">
<summary><strong>Understanding Magnitude Scores</strong></summary>

<p>Each annotation source uses its own criteria to assign a magnitude score
(0&ndash;10) that reflects clinical importance. Higher scores warrant more
attention.</p>

<h3>ClinVar (clinical significance)</h3>
<table style="width:auto; font-size:.85rem; margin:.5rem 0;">
<tr><td><strong>9</strong></td><td>Pathogenic</td></tr>
<tr><td><strong>8</strong></td><td>Pathogenic (single submitter / no assertion criteria)</td></tr>
<tr><td><strong>7</strong></td><td>Likely pathogenic</td></tr>
<tr><td><strong>5</strong></td><td>Uncertain significance / conflicting</td></tr>
<tr><td><strong>4</strong></td><td>Risk factor / drug response / association</td></tr>
<tr><td><strong>3</strong></td><td>Likely benign</td></tr>
<tr><td><strong>1</strong></td><td>Benign</td></tr>
</table>

<h3>PharmGKB (pharmacogenomic evidence)</h3>
<table style="width:auto; font-size:.85rem; margin:.5rem 0;">
<tr><td><strong>9</strong></td><td>Level 1A &mdash; CPIC guideline or FDA label</td></tr>
<tr><td><strong>8</strong></td><td>Level 1B &mdash; strong clinical evidence</td></tr>
<tr><td><strong>6</strong></td><td>Level 2A &mdash; moderate evidence</td></tr>
<tr><td><strong>5</strong></td><td>Level 2B &mdash; moderate (weaker replication)</td></tr>
<tr><td><strong>4</strong></td><td>Level 3 &mdash; low evidence or annotation only</td></tr>
</table>

<h3>GWAS Catalog (trait associations)</h3>
<table style="width:auto; font-size:.85rem; margin:.5rem 0;">
<tr><td><strong>6</strong></td><td>p &lt; 5&times;10<sup>-8</sup> (genome-wide)</td></tr>
<tr><td><strong>4</strong></td><td>p &lt; 5&times;10<sup>-6</sup> (suggestive)</td></tr>
<tr><td><strong>3</strong></td><td>p &lt; 1&times;10<sup>-4</sup> (nominal)</td></tr>
</table>

<h3>SNPedia (community-curated)</h3>
<p>Magnitude is assigned directly by community editors (0&ndash;10 scale).
Higher scores indicate greater clinical or personal relevance as judged by
contributors. See <a href="https://www.snpedia.com/index.php/Magnitude">
SNPedia&rsquo;s magnitude documentation</a> for details.</p>

</details>
"""


_BAD_SIGNIFICANCE = frozenset(
    {
        "clinvar_pathogenic",
        "clinvar_pathogenic/likely_pathogenic",
        "clinvar_likely_pathogenic",
        "clinvar_risk_factor",
        "snpedia_bad",
    }
)

_GOOD_SIGNIFICANCE = frozenset(
    {
        "clinvar_benign",
        "clinvar_benign/likely_benign",
        "clinvar_likely_benign",
        "snpedia_good",
    }
)


def _classify_repute(significance: str) -> str:
    """Derive row border color class from the significance field."""
    sig = significance.lower()
    if sig in _BAD_SIGNIFICANCE:
        return "repute-bad"
    if sig in _GOOD_SIGNIFICANCE:
        return "repute-good"
    return "repute-neutral"


_SORT_SCRIPT = """
<script>
document.addEventListener("DOMContentLoaded",function(){
  var table=document.querySelector(".table-wrap table");
  if(!table)return;
  var headers=table.querySelectorAll("th");
  var tbody=table.querySelector("tbody");
  var dir={};
  headers.forEach(function(th,i){
    th.addEventListener("click",function(){
      var rows=Array.from(tbody.querySelectorAll("tr"));
      dir[i]=dir[i]==="asc"?"desc":"asc";
      rows.sort(function(a,b){
        var ac=a.children[i],bc=b.children[i];
        if(!ac||!bc)return 0;
        var av=ac.getAttribute("data-sort-value")||ac.textContent;
        var bv=bc.getAttribute("data-sort-value")||bc.textContent;
        var an=parseFloat(av),bn=parseFloat(bv);
        if(!isNaN(an)&&!isNaN(bn)){
          return dir[i]==="asc"?an-bn:bn-an;
        }
        return dir[i]==="asc"?av.localeCompare(bv):bv.localeCompare(av);
      });
      rows.forEach(function(r){tbody.appendChild(r);});
      headers.forEach(function(h){
        var arrow=h.querySelector(".sort-arrow");
        if(arrow)arrow.textContent="";
      });
      var arrow=th.querySelector(".sort-arrow");
      if(arrow)arrow.textContent=dir[i]==="asc"?"\\u25B2":"\\u25BC";
    });
  });
});
</script>
"""


def _escape(value: str) -> str:
    return html.escape(value or "", quote=True)


def _hv_nocall_banner(warnings: list[str] | None) -> str:
    if not warnings:
        return ""
    items = "".join(f"<li>{_escape(w)}</li>" for w in warnings)
    return (
        '<div class="notice-warn"><strong>High-value no-calls.</strong> '
        "The following clinically important SNPs returned no genotype call:"
        f"<ul>{items}</ul></div>"
    )


_AM_CSS_CLASS: dict[str, str] = {
    "likely_pathogenic": "am-pathogenic",
    "likely_benign": "am-benign",
    "ambiguous": "am-ambiguous",
}


def _format_am(score: float | None, am_class: str, *, neutral: bool = False) -> str:
    if score is None:
        return "—"
    if neutral:
        return (
            f'<span class="am-score"'
            f' title="AlphaMissense: protein structure impact only">{score:.3f}</span>'
        )
    css = _AM_CSS_CLASS.get(am_class, "am-score")
    return f'<span class="{css}" title="{am_class}">{score:.3f}</span>'


def _format_freq(af: float | None) -> str:
    if af is None:
        return "—"
    pct = af * 100
    if pct < 0.01:
        return "&lt;0.01%"
    return f"{pct:.2f}%"


def _format_cadd(score: float | None) -> str:
    """Format a CADD PHRED score for display."""
    if score is None:
        return "—"
    if score >= 30:
        css, tip = "cadd-high", "top 0.1% most deleterious"
    elif score >= 20:
        css, tip = "cadd-med", "top 1% most deleterious"
    else:
        css, tip = "cadd-low", "top 10% most deleterious" if score >= 10 else ""
    if tip:
        return f'<span class="{css}" title="{tip}">{score:.1f}</span>'
    return f'<span class="{css}">{score:.1f}</span>'


def _row_html(
    a: Annotation,
    css_class: str = "",
    *,
    show_freq: bool = False,
    show_review: bool = True,
    show_am: bool = False,
    show_cadd: bool = False,
) -> str:
    """Render a single annotation as an HTML table row."""
    bar_width = max(1, int(a.magnitude * 8))
    refs_html = ""
    if a.references:
        refs_items = " ".join(_escape(r) for r in a.references)
        refs_html = f'<details class="refs-toggle"><summary>refs</summary>{refs_items}</details>'
    repute = _classify_repute(a.significance)
    classes = [repute]
    if css_class:
        classes.append(css_class)
    tr_open = f'<tr class="{" ".join(classes)}">'
    freq_td = f"<td>{_format_freq(a.allele_frequency)}</td>" if show_freq else ""
    review_td = f"<td>{_escape(a.review_status) or '—'}</td>" if show_review else ""
    am_neutral = a.source == "pharmgkb"
    am_cell = _format_am(a.am_pathogenicity, a.am_class, neutral=am_neutral)
    am_td = f"<td>{am_cell}</td>" if show_am else ""
    cadd_td = f"<td>{_format_cadd(a.cadd_phred)}</td>" if show_cadd else ""
    return (
        f"{tr_open}"
        f'<td class="col-rsid">{_escape(a.rsid)}</td>'
        f'<td class="gene">{_escape(a.gene) or "—"}</td>'
        f'<td><span class="source">{_escape(a.attribution)}</span></td>'
        f"<td>{_escape(a.significance)}</td>"
        f"{review_td}"
        f'<td data-sort-value="{a.magnitude:.1f}">'
        f'<span class="bar" style="width: {bar_width}px;"></span>'
        f"{a.magnitude:.1f}</td>"
        f"<td>{_escape(a.genotype_match)}</td>"
        f"<td>{_escape(a.zygosity)}</td>"
        f"{freq_td}"
        f"{am_td}"
        f"{cadd_td}"
        f'<td class="desc-cell">{_escape(a.condition) or "—"}<br>'
        f'<span class="condition">{_escape(a.description)}</span>{refs_html}</td>'
        "</tr>"
    )


def _removed_row_html(d: dict, *, show_review: bool = True) -> str:
    """Render a removed annotation from a previous report dict."""
    bar_width = max(1, int(d.get("magnitude", 0) * 8))
    mag = d.get("magnitude", 0.0)
    review_td = f"<td>{_escape(d.get('review_status', '')) or '—'}</td>" if show_review else ""
    return (
        '<tr class="diff-removed">'
        f'<td class="col-rsid">{_escape(d.get("rsid", ""))}</td>'
        f'<td class="gene">{_escape(d.get("gene", "")) or "—"}</td>'
        f'<td><span class="source">{_escape(d.get("attribution", ""))}</span></td>'
        f"<td>{_escape(d.get('significance', ''))}</td>"
        f"{review_td}"
        f'<td data-sort-value="{mag:.1f}">'
        f'<span class="bar" style="width: {bar_width}px;"></span>'
        f"{mag:.1f}</td>"
        f"<td>{_escape(d.get('genotype_match', ''))}</td>"
        f'<td class="desc-cell">{_escape(d.get("condition", "")) or "—"}<br>'
        f'<span class="condition">{_escape(d.get("description", ""))}</span></td>'
        "</tr>"
    )


def _license_attributions(annotators_used: list[tuple[str, str | None]]) -> str:
    """Build license attribution HTML from annotator LicenseDescriptors."""
    import logging

    from allelix.annotators import get_annotator_class

    logger = logging.getLogger(__name__)
    parts: list[str] = []
    for name, _version in annotators_used:
        cls = get_annotator_class(name)
        if cls is None:
            logger.warning("No annotator class found for '%s' — attribution omitted", name)
            continue
        desc = cls.license
        source_link = desc.source_url or desc.license_url
        parts.append(
            f" <a href='{html.escape(source_link)}'>"
            f"{html.escape(cls.display_name)}</a>: "
            f"{html.escape(desc.attribution_text)}"
            f" (<a href='{html.escape(desc.license_url)}'>license</a>)"
        )
    return "".join(parts)


def render_html(
    result: AnalysisResult,
    *,
    output_path: Path,
    min_magnitude: float = 0.0,
    category: str | None = None,
    genes: Iterable[str] | None = None,
    source_min_magnitudes: dict[str, float] | None = None,
    title: str = "Allelix Genotype Report",
    diff: DiffResult | None = None,
    high_value_no_calls: list[str] | None = None,
) -> int:
    """Write a self-contained HTML report. Returns the number of annotations rendered."""
    filtered = result.filter(
        min_magnitude=min_magnitude,
        category=category,
        genes=genes,
        source_min_magnitudes=source_min_magnitudes,
    )
    filtered = rollup_gwas_duplicates(filtered)
    annotators_str = ", ".join(
        f"{name} ({version or 'unknown'})" for name, version in result.annotators_used
    )

    build_warn = ""
    diag = result.build_diagnostics
    if diag is not None and diag.mismatch:
        build_warn = (
            '<div class="notice-warn"><strong>Build mismatch.</strong> '
            f"File header claims {_escape(diag.header_build or '')} but position data "
            f"indicates {_escape(diag.detected_build or '')}. "
            f"This report uses {_escape(diag.effective_build)}. "
            "Your provider may have mislabeled the genome build.</div>"
        )

    diff_banner = ""
    new_keys: set[tuple[str, str, str, str]] = set()
    changed_keys: set[tuple[str, str, str, str]] = set()
    if diff is not None:
        from allelix.reports.diff import summarize_diff

        diff_banner = (
            f'<div class="notice"><strong>Diff: </strong>{_escape(summarize_diff(diff))}</div>'
        )
        new_keys = {(a.source, a.rsid, a.condition, a.description) for a in diff.new}
        changed_keys = {
            (c.current.source, c.current.rsid, c.current.condition, c.current.description)
            for c in diff.changed
        }

    floor_note = ""
    if source_min_magnitudes and min_magnitude > 0:
        lower_floors = {
            src: floor for src, floor in source_min_magnitudes.items() if floor < min_magnitude
        }
        if lower_floors:
            parts = ", ".join(
                f"{src.upper()} (mag &ge; {floor:.1f})"
                for src, floor in sorted(lower_floors.items())
            )
            floor_note = (
                '<div class="notice">'
                f"<strong>Source-specific thresholds.</strong> "
                f"The global minimum magnitude is {min_magnitude:.1f}, "
                f"but lower thresholds apply to: {parts}. "
                "Some rows below the global threshold may appear from these sources."
                "</div>"
            )

    has_freq = any(a.allele_frequency is not None for a in filtered)
    show_review = any(a.review_status and a.review_status != "—" for a in filtered)
    has_am = any(a.am_pathogenicity is not None for a in filtered)
    has_cadd = any(a.cadd_phred is not None for a in filtered)

    def _th(label: str) -> str:
        return f'<th>{label}<span class="sort-arrow"></span></th>'

    if filtered or (diff and diff.removed):
        rows_parts: list[str] = []
        for a in filtered:
            key = (a.source, a.rsid, a.condition, a.description)
            if key in new_keys:
                css = "diff-new"
            elif key in changed_keys:
                css = "diff-changed"
            else:
                css = ""
            rows_parts.append(
                _row_html(
                    a,
                    css_class=css,
                    show_freq=has_freq,
                    show_review=show_review,
                    show_am=has_am,
                    show_cadd=has_cadd,
                )
            )
        if diff and diff.removed:
            rows_parts.extend(_removed_row_html(d, show_review=show_review) for d in diff.removed)
        rows_html = "\n".join(rows_parts)
        freq_th = _th("Pop. Freq") if has_freq else ""
        review_th = _th("Review Status") if show_review else ""
        am_th = _th("AM") if has_am else ""
        cadd_th = _th("CADD") if has_cadd else ""
        body = (
            '<div class="table-wrap">'
            "<table>"
            "<thead><tr>"
            f"{_th('rsID')}{_th('Gene')}{_th('Source')}{_th('Significance')}"
            f"{review_th}{_th('Magnitude')}{_th('Genotype')}{_th('Zygosity')}"
            f"{freq_th}"
            f"{am_th}"
            f"{cadd_th}"
            f"{_th('Condition / Description')}"
            "</tr></thead>"
            f"<tbody>{rows_html}</tbody></table>"
            "</div>"
        )
    else:
        body = '<div class="empty">No annotations matched the current filters.</div>'

    summary_cards = "\n".join(
        f'<div class="card"><div class="label">{label}</div><div class="value">{value}</div></div>'
        for label, value in [
            ("Sample", _escape(result.sample_id) or "(unknown)"),
            ("Format", _escape(result.parser_display_name)),
            ("Build", _escape(result.build)),
            ("Variants", f"{result.total_variants:,}"),
            ("Annotations", f"{len(filtered):,}"),
        ]
    )

    document = (
        "<!DOCTYPE html>"
        "<html lang='en'><head><meta charset='utf-8'>"
        f"<title>{_escape(title)}</title>"
        f"<style>{_CSS}</style>"
        "</head><body>"
        f"<h1>{_escape(title)}</h1>"
        f'<p class="subtitle">Source: <code>{_escape(result.file_path.name)}</code> · '
        f"Annotators: {_escape(annotators_str)} · "
        f"Generated {datetime.now(UTC).strftime('%Y-%m-%d %H:%M UTC')}</p>"
        '<div class="notice"><strong>Informational only.</strong> '
        f"{_escape(REGULATORY_NOTICE)}</div>"
        f"{build_warn}"
        f"{_hv_nocall_banner(high_value_no_calls)}"
        f"{_EDUCATION_SECTION}"
        f"{_MAGNITUDE_LEGEND}"
        f"{diff_banner}"
        f'<div class="summary">{summary_cards}</div>'
        f"{floor_note}"
        f"{body}"
        f"{_SORT_SCRIPT}"
        f"<footer>Generated by Allelix v{_escape(__version__)} — "
        "<a href='https://github.com/dial481/allelix'>github.com/dial481/allelix</a>. "
        "All variant classifications attributed to their source databases."
        f"{_license_attributions(result.annotators_used)}</footer>"
        "</body></html>"
    )
    atomic_write_text(output_path, document)
    return len(filtered)
