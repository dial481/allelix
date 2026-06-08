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
.summary {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
  gap: .75rem; margin-bottom: 2rem;
}
.card { background: #f5f5f5; padding: .75rem 1rem; border-radius: 4px; }
.card .label {
  font-size: .8rem; color: #666;
  text-transform: uppercase; letter-spacing: .05em;
}
.card .value { font-size: 1.2rem; font-weight: 600; }
table { width: 100%; border-collapse: collapse; font-size: .9rem; }
th, td {
  text-align: left; padding: .55rem .5rem;
  border-bottom: 1px solid #eee; vertical-align: top;
}
th { background: #fafafa; position: sticky; top: 0; font-weight: 600; }
tr:hover { background: #fcfcfc; }
.rsid { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }
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
.empty { padding: 2rem; text-align: center; color: #666; font-style: italic; }
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


def _format_freq(af: float | None) -> str:
    if af is None:
        return "—"
    pct = af * 100
    if pct < 0.01:
        return "&lt;0.01%"
    return f"{pct:.2f}%"


def _row_html(a: Annotation, css_class: str = "", *, show_freq: bool = False) -> str:
    bar_width = max(1, int(a.magnitude * 8))
    refs_html = ""
    if a.references:
        refs_html = (
            '<div class="condition">refs: ' + " ".join(_escape(r) for r in a.references) + "</div>"
        )
    tr_open = f'<tr class="{css_class}">' if css_class else "<tr>"
    freq_td = f"<td>{_format_freq(a.allele_frequency)}</td>" if show_freq else ""
    return (
        f"{tr_open}"
        f'<td class="rsid">{_escape(a.rsid)}</td>'
        f'<td class="gene">{_escape(a.gene) or "—"}</td>'
        f'<td><span class="source">{_escape(a.attribution)}</span></td>'
        f"<td>{_escape(a.significance)}</td>"
        f"<td>{_escape(a.review_status) or '—'}</td>"
        f'<td><span class="bar" style="width: {bar_width}px;"></span>'
        f"{a.magnitude:.1f}</td>"
        f"<td>{_escape(a.genotype_match)}</td>"
        f"{freq_td}"
        f"<td>{_escape(a.condition) or '—'}<br>"
        f'<span class="condition">{_escape(a.description)}</span>{refs_html}</td>'
        "</tr>"
    )


def _removed_row_html(d: dict) -> str:
    """Render a removed annotation from a previous report dict."""
    bar_width = max(1, int(d.get("magnitude", 0) * 8))
    return (
        '<tr class="diff-removed">'
        f'<td class="rsid">{_escape(d.get("rsid", ""))}</td>'
        f'<td class="gene">{_escape(d.get("gene", "")) or "—"}</td>'
        f'<td><span class="source">{_escape(d.get("attribution", ""))}</span></td>'
        f"<td>{_escape(d.get('significance', ''))}</td>"
        f"<td>{_escape(d.get('review_status', '')) or '—'}</td>"
        f'<td><span class="bar" style="width: {bar_width}px;"></span>'
        f"{d.get('magnitude', 0.0):.1f}</td>"
        f"<td>{_escape(d.get('genotype_match', ''))}</td>"
        f"<td>{_escape(d.get('condition', '')) or '—'}<br>"
        f'<span class="condition">{_escape(d.get("description", ""))}</span></td>'
        "</tr>"
    )


_LICENSE_ATTRIBUTIONS: dict[str, str] = {
    "pharmgkb": (
        " Pharmacogenomic annotations sourced from"
        " <a href='https://www.pharmgkb.org'>PharmGKB</a>,"
        " used under <a href='https://creativecommons.org/licenses/by-sa/4.0/'>CC BY-SA 4.0</a>."
    ),
    "snpedia": (
        " SNPedia annotations sourced from"
        " <a href='https://www.snpedia.com'>SNPedia</a>,"
        " used under"
        " <a href='https://creativecommons.org/licenses/by-nc-sa/3.0/us/'>CC BY-NC-SA 3.0 US</a>."
    ),
    "gnomad": (
        " Population frequencies sourced from"
        " <a href='https://gnomad.broadinstitute.org'>gnomAD</a>,"
        " used under <a href='https://opendatacommons.org/licenses/odbl/1-0/'>ODbL v1.0</a>."
    ),
}


def _license_attributions(annotators_used: list[tuple[str, str | None]]) -> str:
    """Build license attribution HTML for annotators that require it."""
    names = {name for name, _version in annotators_used}
    parts = [text for key, text in _LICENSE_ATTRIBUTIONS.items() if key in names]
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

    has_freq = any(a.allele_frequency is not None for a in filtered)

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
            rows_parts.append(_row_html(a, css_class=css, show_freq=has_freq))
        if diff and diff.removed:
            rows_parts.extend(_removed_row_html(d) for d in diff.removed)
        rows_html = "\n".join(rows_parts)
        freq_th = "<th>Pop. Freq</th>" if has_freq else ""
        body = (
            "<table>"
            "<thead><tr>"
            "<th>rsID</th><th>Gene</th><th>Source</th><th>Significance</th>"
            f"<th>Review Status</th><th>Magnitude</th><th>Genotype</th>"
            f"{freq_th}"
            "<th>Condition / Description</th>"
            "</tr></thead>"
            f"<tbody>{rows_html}</tbody></table>"
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
        f"{diff_banner}"
        f'<div class="summary">{summary_cards}</div>'
        f"{body}"
        f"<footer>Generated by Allelix v{_escape(__version__)} — "
        "<a href='https://github.com/dial481/allelix'>github.com/dial481/allelix</a>. "
        "All variant classifications attributed to their source databases."
        f"{_license_attributions(result.annotators_used)}</footer>"
        "</body></html>"
    )
    atomic_write_text(output_path, document)
    return len(filtered)
