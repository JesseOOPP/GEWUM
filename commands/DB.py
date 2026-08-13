"""GEWUM DB command (P1: read-only inspection + batch CIF export).

Sub-commands:
  info     -- per-DB metadata: row counts by stage, energy range, schema
  list     -- table/csv/json view of filtered rows
  stats    -- GROUP BY pivot (1 or 2 dimensions: formula / stage / sg)
  show     -- dump cif_content for a single (sg, name) row
  export   -- batch export filtered CIFs to a directory or a zip archive

This command never writes to any structures.db: all sqlite connections
are opened mode=ro (read-only URI).
"""

from __future__ import annotations

import os
import sys

from ..src.common.db_admin import query as dbq
from ..src.common.db_admin import formatter as dbf
from ..src.common.db_admin import stats as dbs
from ..src.common.db_admin import exporter as dbe

def setup_args(parser):
    parser.description = (
        "GEWUM DB - read-only inspection and batch CIF export from "
        "structures.db files."
    )
    sub = parser.add_subparsers(dest="db_cmd", metavar="<sub-command>")
    sub.required = False  

    _add_common_source_args(sub.add_parser("info",
        help="Show per-DB metadata: schema, row counts by stage, energy range."))

    p_list = sub.add_parser("list",
        help="List rows matching a filter (table/csv/json/tsv).")
    _add_common_source_args(p_list)
    _add_common_filter_args(p_list)
    p_list.add_argument("--format", default="table",
                        choices=["table", "csv", "tsv", "json"],
                        help="Output format (default: table).")
    p_list.add_argument("--columns", default=None,
                        help="Comma-separated columns to display "
                             "(aliases epa/sg/name allowed).")

    p_stats = sub.add_parser("stats",
        help="Aggregate counts (and epa min/avg/max) by 1 or 2 dimensions.")
    _add_common_source_args(p_stats)
    _add_common_filter_args(p_stats)
    p_stats.add_argument("--by", required=True,
                         help="Dimension(s) to group by, comma-separated. "
                              "Allowed: formula, stage, sg. Max 2.")
    p_stats.add_argument("--no-epa", action="store_true",
                         help="Suppress epa min/avg/max columns (1-D mode).")

    p_show = sub.add_parser("show",
        help="Dump cif_content for a single (sg_number, cif_name) row.")
    _add_common_source_args(p_show, allow_root=False)  # show targets one DB
    p_show.add_argument("--sg", type=int, required=True, help="Space group number.")
    p_show.add_argument("--name", required=True, help="cif_name (exact).")
    p_show.add_argument("--initial", action="store_true",
                        help="Dump the pre-relax CIF snapshot (cif_content_initial); "
                             "falls back to live cif_content for rows that have "
                             "never been relaxed.")
    p_show.add_argument("--to", default=None,
                        help="Write CIF to this path; default: stdout.")

    p_exp = sub.add_parser("export",
        help="Batch export filtered CIFs to a directory or zip archive.")
    _add_common_source_args(p_exp)
    _add_common_filter_args(p_exp)
    out_grp = p_exp.add_mutually_exclusive_group(required=True)
    out_grp.add_argument("--out", default=None,
                         help="Output directory for exported .cif files.")
    out_grp.add_argument("--zip", dest="zip_path", default=None,
                         help="Output zip archive path.")
    p_exp.add_argument("--layout", default="by-formula",
                       choices=list(dbe.LAYOUTS),
                       help="Sub-directory layout (default: by-formula).")
    p_exp.add_argument("--name-template",
                       default="{name}.cif",
                       help="Filename template; placeholders: "
                            "{formula} {sg} {name} {cif} {stage} {epa} {energy}. "
                            "Default: '{name}.cif'.")
    p_exp.add_argument("--workers", type=int, default=1,
                       help="Parallel writer threads (dir mode only).")
    p_exp.add_argument("--overwrite", default="skip",
                       choices=list(dbe.OVERWRITE_MODES),
                       help="Behaviour when target file already exists.")
    p_exp.add_argument("--manifest", default=None,
                       help="Optional CSV manifest path.")
    p_exp.add_argument("--max-files", type=int, default=5000,
                       help="Safety cap on number of files exported (default 5000).")
    p_exp.add_argument("--confirm", action="store_true",
                       help="Required when match count > --max-files.")
    p_exp.add_argument("--dry-run", action="store_true",
                       help="Print what would be written; do not write.")
    p_exp.add_argument("--source", default="relaxed",
                       choices=["relaxed", "initial"],
                       help="Which CIF text to export per row: 'relaxed' (default, "
                            "current cif_content) or 'initial' (pre-relax snapshot).")


def _add_common_source_args(parser, allow_root: bool = True):
    parser.add_argument("--db", action="append", default=None,
                        help="Path to a structures.db file. May be repeated.")
    if allow_root:
        parser.add_argument("--root", action="append", default=None,
                            help="Directory; auto-discovers <root>/<formula>/structures.db. "
                                 "May be repeated.")


def _add_common_filter_args(parser):
    parser.add_argument("--formula", default=None,
                        help="Glob applied to formula folder name (e.g. 'Na*Cl*').")
    parser.add_argument("--sg", default=None,
                        help="Space-group filter: '225', '225,227', '100-150' or mixed.")
    parser.add_argument("--name", default=None,
                        help="Glob applied to cif_name (e.g. 'xtal_*').")
    parser.add_argument("--stage", default=None,
                        help="Comma-separated stages (e.g. 'relaxed,bond_mis').")
    parser.add_argument("--energy-min", type=float, default=None)
    parser.add_argument("--energy-max", type=float, default=None)
    parser.add_argument("--epa-min", type=float, default=None,
                        help="Min energy_per_atom.")
    parser.add_argument("--epa-max", type=float, default=None,
                        help="Max energy_per_atom.")
    parser.add_argument("--order-by", default=None,
                        help="Order rows by 'col[:asc|desc]' (aliases ok).")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--offset", type=int, default=0)

def execute(args, remaining_args=None):
    cmd = getattr(args, "db_cmd", None)
    if not cmd:
        print("[ERROR] gewum DB requires a sub-command "
              "(info | list | stats | show | export).", file=sys.stderr)
        print("Run 'gewum DB -h' for help.", file=sys.stderr)
        sys.exit(2)

    try:
        if cmd == "info":
            _cmd_info(args)
        elif cmd == "list":
            _cmd_list(args)
        elif cmd == "stats":
            _cmd_stats(args)
        elif cmd == "show":
            _cmd_show(args)
        elif cmd == "export":
            _cmd_export(args)
        else:
            print(f"[ERROR] unknown DB sub-command: {cmd}", file=sys.stderr)
            sys.exit(2)
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

def _resolve_dbs(args):
    explicit = getattr(args, "db", None) or []
    roots = getattr(args, "root", None) or []
    if not explicit and not roots:
        raise ValueError(
            "no DB source supplied: pass --db <structures.db> (repeatable) "
            "and/or --root <dir> (repeatable)."
        )
    dbs = dbq.discover_dbs(explicit_dbs=explicit, roots=roots)
    if not dbs:
        raise ValueError("no structures.db found from the given sources.")
    return dbs


def _build_spec(args) -> dbq.FilterSpec:
    stages = None
    if getattr(args, "stage", None):
        stages = [s.strip() for s in args.stage.split(",") if s.strip()]
    return dbq.FilterSpec(
        formula_glob=getattr(args, "formula", None),
        sg_spec=getattr(args, "sg", None),
        name_glob=getattr(args, "name", None),
        stages=stages,
        energy_min=getattr(args, "energy_min", None),
        energy_max=getattr(args, "energy_max", None),
        epa_min=getattr(args, "epa_min", None),
        epa_max=getattr(args, "epa_max", None),
        order_by=getattr(args, "order_by", None),
        limit=getattr(args, "limit", None),
        offset=getattr(args, "offset", 0) or 0,
    )

def _cmd_info(args):
    dbs = _resolve_dbs(args)
    print(f"Discovered {len(dbs)} structures.db file(s):\n")
    for formula_name, db_path in dbs:
        size = os.path.getsize(db_path)
        conn = dbq.open_ro(db_path)
        try:
            cols = dbq.list_columns(conn)
            total = conn.execute("SELECT COUNT(*) FROM structures").fetchone()[0]
            stage_rows = conn.execute(
                "SELECT stage, COUNT(*) FROM structures GROUP BY stage ORDER BY stage"
            ).fetchall()
            sg_lo, sg_hi, sg_distinct = conn.execute(
                "SELECT MIN(sg_number), MAX(sg_number), COUNT(DISTINCT sg_number) "
                "FROM structures"
            ).fetchone()
            epa_lo, epa_hi = conn.execute(
                "SELECT MIN(energy_per_atom), MAX(energy_per_atom) "
                "FROM structures WHERE stage='relaxed'"
            ).fetchone()
            initial_missing = None
            if "cif_content_initial" in cols:
                initial_missing = conn.execute(
                    "SELECT COUNT(*) FROM structures "
                    "WHERE stage='relaxed' AND cif_content_initial IS NULL"
                ).fetchone()[0]
        finally:
            conn.close()

        print(f"== {formula_name} ==")
        print(f"  path     : {db_path}")
        print(f"  size     : {size/1024:.1f} KB")
        print(f"  columns  : {', '.join(cols)}")
        print(f"  rows     : {total}")
        if stage_rows:
            parts = [f"{r['stage'] if isinstance(r, dict) else r[0]}={r[1]}"
                     for r in stage_rows]
            print(f"  by stage : {'  '.join(parts)}")
        if sg_lo is not None:
            print(f"  sg range : {sg_lo} .. {sg_hi}  ({sg_distinct} distinct)")
        if epa_lo is not None:
            print(f"  epa(relaxed): min={epa_lo:.6g}  max={epa_hi:.6g}")
        if "cif_content_initial" not in cols:
            print("  initial CIF: column missing -- legacy DB, pre-relax "
                  "snapshots unavailable for any row")
        elif initial_missing:
            print(f"  initial CIF: {initial_missing} relaxed row(s) lack "
                  f"the pre-relax snapshot (relaxed before snapshot column "
                  f"was introduced; data unrecoverable)")
        else:
            print("  initial CIF: snapshot intact for all relaxed rows")
        print()


def _cmd_list(args):
    dbs = _resolve_dbs(args)
    spec = _build_spec(args)
    columns = None
    if args.columns:
        columns = [c.strip() for c in args.columns.split(",") if c.strip()]
    rows = dbq.iter_rows(dbs, spec, columns=columns)
    dbf.render(rows, fmt=args.format, columns=columns)


def _cmd_stats(args):
    db_list = _resolve_dbs(args)
    spec = _build_spec(args)
    by = [b.strip() for b in args.by.split(",") if b.strip()]
    keys, agg = dbs.aggregate(db_list, spec, by)
    show_epa = (len(by) == 1) and (not args.no_epa)
    print(dbs.render_stats(by, keys, agg, show_epa=show_epa))


def _cmd_show(args):
    explicit = getattr(args, "db", None) or []
    if not explicit:
        raise ValueError("'show' requires --db <structures.db>")
    dbs_list = dbq.discover_dbs(explicit_dbs=explicit, roots=None)
    if not dbs_list:
        raise ValueError("no DB found from --db argument")
    if len(dbs_list) > 1:
        raise ValueError("'show' accepts exactly one --db; "
                         f"got {len(dbs_list)}")
    formula, db_path = dbs_list[0]
    conn = dbq.open_ro(db_path)
    try:
        cols = dbq.list_columns(conn)
        has_snapshot = "cif_content_initial" in cols
        if args.initial:
            if has_snapshot:
                expr = "COALESCE(cif_content_initial, cif_content)"
            else:
                expr = "cif_content"
            sql = f"SELECT {expr} AS cif, stage FROM structures " \
                  "WHERE sg_number=? AND cif_name=?"
        else:
            sql = "SELECT cif_content AS cif, stage FROM structures " \
                  "WHERE sg_number=? AND cif_name=?"
        row = conn.execute(sql, (args.sg, args.name)).fetchone()
    finally:
        conn.close()
    if not row:
        raise ValueError(
            f"no row in {db_path} with sg_number={args.sg} cif_name={args.name!r}"
        )
    content = row["cif"] or ""
    if args.initial and row["stage"] == "relaxed" and not has_snapshot:
        print("[show] WARNING: --initial requested on a legacy DB without the "
              "snapshot column; this row is already relaxed so the dumped CIF "
              "is the RELAXED structure, not the initial one.", file=sys.stderr)
    if args.to:
        os.makedirs(os.path.dirname(os.path.abspath(args.to)) or ".", exist_ok=True)
        with open(args.to, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"[show] wrote {args.to}  ({len(content)} chars, "
              f"source={'initial' if args.initial else 'relaxed'})")
    else:
        sys.stdout.write(content)
        if not content.endswith("\n"):
            sys.stdout.write("\n")


def _cmd_export(args):
    db_list = _resolve_dbs(args)
    spec = _build_spec(args)
    summary = dbe.export(
        dbs=db_list,
        spec=spec,
        out_dir=args.out,
        zip_path=args.zip_path,
        layout=args.layout,
        name_template=args.name_template,
        workers=max(1, int(args.workers or 1)),
        overwrite=args.overwrite,
        manifest=args.manifest,
        max_files=int(args.max_files),
        confirm_above_max=bool(args.confirm),
        dry_run=bool(args.dry_run),
        source=args.source,
    )
    if summary.get("failed"):
        sys.exit(1)
