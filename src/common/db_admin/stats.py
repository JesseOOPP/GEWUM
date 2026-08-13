"""Aggregation (GROUP BY pivot) for GEWUM DB admin.

Supports one or two grouping dimensions (--by formula,stage etc.) and
the standard aggregates needed for stats: count, plus min/avg/max of
energy_per_atom for the relaxed subset.

Aggregation is performed per DB then merged in Python to keep the SQL
trivial and avoid attaching cross-DB databases.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

from .query import (
    FilterSpec, list_columns, open_ro,
)

DIMENSION_SQL = {
    "formula": "<<formula_literal>>",  
    "stage":   "stage",
    "sg":      "sg_number",
}


def _check_dim(dim: str) -> None:
    if dim not in DIMENSION_SQL:
        raise ValueError(
            f"unsupported --by dimension '{dim}' "
            f"(allowed: {', '.join(DIMENSION_SQL)})"
        )


def aggregate(
    dbs: Sequence[Tuple[str, str]],
    spec: FilterSpec,
    by: Sequence[str],
) -> Tuple[List[Tuple], Dict[Tuple, dict]]:
    """Group rows by the given dimensions and aggregate basic stats.

    Args:
        dbs:  list of (formula_name, db_path).
        spec: pre-aggregation row filter (LIMIT/ORDER are ignored here).
        by:   1 or 2 dimension names; allowed values are keys of DIMENSION_SQL.

    Returns:
        (key_order, agg_map)
          key_order: sorted list of group keys (1- or 2-tuples) in display order
          agg_map:   dict[group_key] -> {'count', 'epa_min', 'epa_avg', 'epa_max'}
    """
    if not by:
        raise ValueError("--by must specify at least one dimension")
    if len(by) > 2:
        raise ValueError("--by accepts at most 2 dimensions")
    for d in by:
        _check_dim(d)

    agg: Dict[Tuple, dict] = defaultdict(
        lambda: {"count": 0, "epa_sum": 0.0, "epa_cnt": 0,
                 "epa_min": None, "epa_max": None}
    )

    pure_spec = FilterSpec(
        formula_glob=spec.formula_glob,
        sg_spec=spec.sg_spec,
        name_glob=spec.name_glob,
        stages=spec.stages,
        energy_min=spec.energy_min,
        energy_max=spec.energy_max,
        epa_min=spec.epa_min,
        epa_max=spec.epa_max,
    )

    sql_dims = [d for d in by if d != "formula"]
    sql_select = [DIMENSION_SQL[d] for d in sql_dims] if sql_dims else []

    for formula_name, db_path in dbs:
        if not pure_spec.matches_formula(formula_name):
            continue
        conn = open_ro(db_path)
        try:
            available = list_columns(conn)
            where, params = pure_spec.to_sql(available)

            agg_select = (
                "COUNT(*), MIN(energy_per_atom), MAX(energy_per_atom), "
                "SUM(energy_per_atom), "
                "SUM(CASE WHEN energy_per_atom IS NOT NULL THEN 1 ELSE 0 END)"
            )
            if sql_dims:
                sql = (
                    f"SELECT {', '.join(sql_select)}, {agg_select} "
                    f"FROM structures {where} "
                    f"GROUP BY {', '.join(sql_select)}"
                )
            else:
                sql = f"SELECT {agg_select} FROM structures {where}"

            for row in conn.execute(sql, params):
                sql_idx = 0
                key_vals = []
                for d in by:
                    if d == "formula":
                        key_vals.append(formula_name)
                    else:
                        key_vals.append(row[sql_idx])
                        sql_idx += 1
                base = sql_idx  
                cnt = row[base]
                if cnt == 0:
                    continue
                vmin = row[base + 1]
                vmax = row[base + 2]
                vsum = row[base + 3] or 0.0
                vcnt = row[base + 4] or 0
                bucket = agg[tuple(key_vals)]
                bucket["count"] += cnt
                if vmin is not None:
                    bucket["epa_min"] = vmin if bucket["epa_min"] is None else min(bucket["epa_min"], vmin)
                if vmax is not None:
                    bucket["epa_max"] = vmax if bucket["epa_max"] is None else max(bucket["epa_max"], vmax)
                bucket["epa_sum"] += vsum
                bucket["epa_cnt"] += vcnt
        finally:
            conn.close()

    for k, b in agg.items():
        b["epa_avg"] = (b["epa_sum"] / b["epa_cnt"]) if b["epa_cnt"] else None
        b.pop("epa_sum", None)
        b.pop("epa_cnt", None)

    keys = sorted(agg.keys(), key=lambda k: tuple("" if x is None else x for x in k))
    return keys, dict(agg)

def render_stats(
    by: Sequence[str],
    keys: List[Tuple],
    agg: Dict[Tuple, dict],
    show_epa: bool = True,
) -> str:
    if not keys:
        return "(no rows match the filter)"

    if len(by) == 1:
        headers = [by[0], "count"]
        if show_epa:
            headers += ["epa_min", "epa_avg", "epa_max"]
        rows = []
        total = 0
        for k in keys:
            d = agg[k]
            row = [str(k[0]), str(d["count"])]
            if show_epa:
                row += [
                    "" if d["epa_min"] is None else f"{d['epa_min']:.6g}",
                    "" if d["epa_avg"] is None else f"{d['epa_avg']:.6g}",
                    "" if d["epa_max"] is None else f"{d['epa_max']:.6g}",
                ]
            rows.append(row)
            total += d["count"]
        rows.append(["TOTAL", str(total)] + [""] * (len(headers) - 2))
        return _table(headers, rows)

    row_keys = sorted({k[0] for k in keys}, key=lambda x: ("" if x is None else x))
    col_keys = sorted({k[1] for k in keys}, key=lambda x: ("" if x is None else x))
    headers = [by[0]] + [str(c) for c in col_keys] + ["TOTAL"]
    rows = []
    col_totals = [0] * len(col_keys)
    grand_total = 0
    for r in row_keys:
        line = [str(r)]
        rt = 0
        for j, c in enumerate(col_keys):
            v = agg.get((r, c), {}).get("count", 0)
            line.append(str(v))
            rt += v
            col_totals[j] += v
        line.append(str(rt))
        rows.append(line)
        grand_total += rt
    rows.append(["TOTAL"] + [str(t) for t in col_totals] + [str(grand_total)])
    return _table(headers, rows)


def _table(headers: List[str], rows: List[List[str]]) -> str:
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            if len(cell) > widths[i]:
                widths[i] = len(cell)

    def line(vals):
        return "  ".join(v.ljust(widths[i]) for i, v in enumerate(vals))

    out = [line(headers), line(["-" * w for w in widths])]
    out.extend(line(r) for r in rows)
    return "\n".join(out)
