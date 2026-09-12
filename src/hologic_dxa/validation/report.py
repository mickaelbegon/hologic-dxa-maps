"""HTML validation report generator."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from hologic_dxa.validation.mass_conservation import ValidationReport


_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>DXA Validation Report</title>
  <style>
    body {{ font-family: system-ui, sans-serif; margin: 2rem; color: #222; }}
    h1 {{ color: #1a3c5e; }}
    .pass {{ color: #2e7d32; font-weight: bold; }}
    .fail {{ color: #c62828; font-weight: bold; }}
    .warn {{ color: #e65100; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; }}
    th, td {{ border: 1px solid #ccc; padding: 0.5rem 1rem; text-align: left; }}
    th {{ background: #1a3c5e; color: white; }}
    tr:nth-child(even) {{ background: #f5f5f5; }}
    .meta {{ font-size: 0.9rem; color: #555; margin-top: 2rem; }}
    .disclaimer {{ background: #fff8e1; border-left: 4px solid #f9a825;
                   padding: 1rem; margin: 1rem 0; }}
  </style>
</head>
<body>
  <h1>DXA Quantitative Map Validation Report</h1>
  <p>Generated: {timestamp}</p>
  <p>Pipeline version: {pipeline_version}</p>

  <div class="disclaimer">
    <strong>Important:</strong> Acceptance thresholds used in this report are
    <em>provisional</em> and have not been validated against phantom data.
    A "PASS" result is conditional on the quality of the input calibration.
    Do not interpret passing results as metrological certification.
  </div>

  <h2>Overall result: <span class="{overall_class}">{overall_text}</span></h2>

  <h2>Regional comparisons</h2>
  <table>
    <thead>
      <tr>
        <th>Region</th>
        <th>Integrated [g]</th>
        <th>Reference [g]</th>
        <th>Absolute error [g]</th>
        <th>Relative error</th>
        <th>N valid pixels</th>
        <th>Result</th>
      </tr>
    </thead>
    <tbody>
{comparison_rows}
    </tbody>
  </table>

  {warnings_section}

  <div class="meta">
    <h3>Thresholds used</h3>
    <pre>{thresholds}</pre>
    <h3>Provenance</h3>
    <pre>{provenance}</pre>
  </div>
</body>
</html>
"""


def generate_html_report(report: "ValidationReport", output_path: Path) -> None:
    """Write an HTML validation report to output_path."""
    import json

    rows = []
    for comp in report.comparisons:
        result_class = "pass" if comp.passed else "fail"
        result_text = "PASS" if comp.passed else "FAIL"
        rows.append(
            f"      <tr>"
            f"<td>{comp.region_name}</td>"
            f"<td>{comp.integrated_g:.1f}</td>"
            f"<td>{comp.reference_g:.1f}</td>"
            f"<td>{comp.absolute_error_g:.1f}</td>"
            f"<td>{comp.relative_error:+.2%}</td>"
            f"<td>{comp.n_valid_pixels}</td>"
            f"<td class='{result_class}'>{result_text}</td>"
            f"</tr>"
        )

    warnings_html = ""
    if report.warnings:
        items = "".join(f"<li class='warn'>{w}</li>" for w in report.warnings)
        warnings_html = f"<h2>Warnings</h2><ul>{items}</ul>"

    html = _HTML_TEMPLATE.format(
        timestamp=report.timestamp_utc.isoformat(),
        pipeline_version=report.pipeline_version,
        overall_class="pass" if report.overall_passed else "fail",
        overall_text="PASS" if report.overall_passed else "FAIL",
        comparison_rows="\n".join(rows) if rows else "      <tr><td colspan='7'>No comparisons available.</td></tr>",
        warnings_section=warnings_html,
        thresholds=json.dumps(report.thresholds_used, indent=2),
        provenance=json.dumps(report.provenance, indent=2, default=str),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
