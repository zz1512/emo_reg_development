"""Small shared IO helpers for read-only reviewer-table extractors."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.fig03.gam_pygam import result_tables

ROOT = Path(__file__).resolve().parents[2]


def project_path(path):
    path = Path(path)
    return (path if path.is_absolute() else ROOT / path).resolve()


def read_csv(path):
    return pd.read_csv(path, encoding="utf-8-sig", float_precision="round_trip",
                       dtype={"name": str, "sub_id": str, "site_id": str, "outer_site": str})


def required(frame, columns, source):
    if missing := set(columns) - set(frame):
        raise ValueError(f"{source}: missing columns {sorted(missing)}")


def boolean(value):
    if pd.isna(value):
        return None
    text = str(value).strip().lower()
    if text in ("true", "1", "1.0"):
        return True
    if text in ("false", "0", "0.0"):
        return False
    raise ValueError(f"Invalid boolean value: {value!r}")


def close(actual, expected, description):
    if not np.allclose(np.asarray(actual, dtype=float), np.asarray(expected, dtype=float),
                       rtol=1e-8, atol=1e-10, equal_nan=True):
        raise ValueError(f"Inconsistent saved results: {description}")


def number(value):
    if pd.isna(value):
        return "NA"
    if np.isinf(value):
        return "overflow; see log BF"
    if value != 0 and (abs(value) < .001 or abs(value) >= 10000):
        return f"{value:.2e}"
    return f"{value:.3f}"


def interval(low, high, integer=False):
    if pd.isna(low) or pd.isna(high):
        return "NA"
    render = (lambda value: str(int(value))) if integer else number
    return render(low) if low == high else f"{render(low)} to {render(high)}"


def markdown_table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers),
                      *("| " + " | ".join(map(str, row)) + " |" for row in rows)])


def export(source, output, kind, frames, text, sections, sources, *, mock=False, docx=True, extra=None):
    """Write only an owned report directory; source analysis files stay untouched."""
    source, output = source.resolve(), output.resolve()
    if output == source or output in source.parents:
        raise ValueError("The report directory must be separate from the source analysis")
    manifest_path = output / "summary_metadata.json"
    previous = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    if output.exists() and any(output.iterdir()) and previous.get("analysis_kind") != kind:
        raise ValueError("Nonempty output is not a directory owned by this report extractor")
    document = None
    if docx:
        document = result_tables._document(sections[0], sections[1], mock)
        for title, headers, rows, widths, note in sections[2]:
            document.add_paragraph(title).runs[0].bold = True
            table = result_tables._three_line_table(document, headers, widths)
            for row in rows:
                result_tables._append_row(table, row)
            document.add_paragraph(note)
    output.mkdir(parents=True, exist_ok=True)
    for name, frame in frames.items():
        frame.to_csv(output / name, index=False, encoding="utf-8-sig", float_format="%.17g")
    if mock:
        text = "**SIMULATED DATA ONLY：仅为模拟，不可填写论文。**\n\n" + text
    (output / "result_summary.md").write_text(text + "\n", encoding="utf-8")
    if document is not None:
        document.save(output / "response_tables.docx")
    elif "response_tables.docx" in previous.get("outputs", []) and (output / "response_tables.docx").is_file():
        (output / "response_tables.docx").unlink()
    result = dict(analysis_kind=kind, simulated=mock,
                  source_sha256={str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                  row_counts={name: len(frame) for name, frame in frames.items()},
                  outputs=[*frames, "result_summary.md", "summary_metadata.json"] + (["response_tables.docx"] if docx else []),
                  source_changes="none; no model refitting, reintegration, resampling or multiple-comparison correction")
    result.update(extra or {})
    manifest_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Report: {output}", flush=True)
    return result
