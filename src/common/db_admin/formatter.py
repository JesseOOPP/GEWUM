"""Output formatters for GEWUM DB admin query results."""

from __future__ import annotations

import csv
import io
import json
import sys
from typing import Iterable, List, Optional, Sequence

DEFAULT_COLUMNS = (
    "formula_name", "sg_number", "cif_name", "stage",
    "energy_per_atom", "energy", "final_pressure", "enthalpy_per_atom",
)


def _format_value(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        if abs(v) >= 1000 or (v != 0 and abs(v) < 1e-3):
            return f"{v:.4e}"
        return f"{v:.6g}"
    if isinstance(v, (bytes, bytearray)):
        return f"<bytes len={len(v)}>"
    return str(v)


def _resolve_columns(rows: List[dict], requested: Optional[Sequence[str]]) -> List[str]:
    if requested:
        return list(requested)
    if not rows:
        return list(DEFAULT_COLUMNS)
    seen = set()
    cols: List[str] = []
    for c in DEFAULT_COLUMNS:
        if any(c in r for r in rows) and c not in seen:
            cols.append(c)
            seen.add(c)
    extras = set()
    for r in rows:
        extras.update(r.keys())
    for c in sorted(extras):
        if c in seen:
            continue
        if c in ("cif_content", "db_path", "created_at"):
            continue
        cols.append(c)
        seen.add(c)
    return cols

def _render_table(rows: List[dict], cols: List[str]) -> str:
    if not rows:
        return "(no rows)"
    string_rows = [[_format_value(r.get(c)) for c in cols] for r in rows]
    widths = [len(c) for c in cols]
    for row in string_rows:
        for i, cell in enumerate(row):
            if len(cell) > widths[i]:
                widths[i] = len(cell)

    def fmt_line(values):
        return "  ".join(v.ljust(widths[i]) for i, v in enumerate(values))

    out = [fmt_line(cols), fmt_line(["-" * w for w in widths])]
    out.extend(fmt_line(r) for r in string_rows)
    return "\n".join(out)


def _render_csv(rows: List[dict], cols: List[str], delimiter: str = ",") -> str:
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter=delimiter, lineterminator="\n")
    writer.writerow(cols)
    for r in rows:
        writer.writerow([_format_value(r.get(c)) for c in cols])
    return buf.getvalue()


def _render_json(rows: List[dict], cols: List[str]) -> str:
    out = [{c: r.get(c) for c in cols} for r in rows]
    return json.dumps(out, indent=2, default=str, ensure_ascii=False)


def render(
    rows: Iterable[dict],
    fmt: str = "table",
    columns: Optional[Sequence[str]] = None,
    out=None,
) -> None:
    """Render query results to a stream (stdout by default)."""
    rows = list(rows)
    cols = _resolve_columns(rows, columns)
    if fmt == "table":
        text = _render_table(rows, cols)
    elif fmt == "csv":
        text = _render_csv(rows, cols, ",")
    elif fmt == "tsv":
        text = _render_csv(rows, cols, "\t")
    elif fmt == "json":
        text = _render_json(rows, cols)
    else:
        raise ValueError(f"Unknown format: {fmt!r} (choose table/csv/tsv/json)")

    stream = out if out is not None else sys.stdout
    stream.write(text)
    if not text.endswith("\n"):
        stream.write("\n")
