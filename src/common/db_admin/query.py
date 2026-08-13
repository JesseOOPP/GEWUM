"""Read-only query layer for GEWUM DB admin.

Builds parameterised SQL from a FilterSpec and yields rows from one or
more structures.db files. All connections are opened in read-only mode
(URI form 'file:...?mode=ro'); inserting any write statement would raise
sqlite3.OperationalError on the connection.
"""

from __future__ import annotations

import os
import sqlite3
import fnmatch
from dataclasses import dataclass, field
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple

DB_NAME = "structures.db"

CORE_COLUMNS = (
    "sg_number", "cif_name", "stage", "energy", "energy_per_atom",
    "formula", "created_at",
)
OPTIONAL_COLUMNS = (
    "final_pressure", "enthalpy_per_atom", "corrected_enthalpy_per_atom",
)
COLUMN_ALIASES = {
    "epa": "energy_per_atom",
    "name": "cif_name",
    "sg":   "sg_number",
    "p":    "final_pressure",
    "h":    "enthalpy_per_atom",
}


def _resolve_column(name: str) -> str:
    return COLUMN_ALIASES.get(name, name)

def open_ro(db_path: str) -> sqlite3.Connection:
    """Open a structures.db read-only.

    Using the URI form with mode=ro guarantees that any write attempt
    raises sqlite3.OperationalError ('attempt to write a readonly database').
    """
    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"structures.db not found: {db_path}")
    uri = "file:" + os.path.abspath(db_path).replace("\\", "/") + "?mode=ro"
    conn = sqlite3.connect(uri, uri=True)
    conn.row_factory = sqlite3.Row
    return conn


def list_columns(conn: sqlite3.Connection) -> List[str]:
    return [r["name"] for r in conn.execute("PRAGMA table_info(structures)").fetchall()]

def discover_dbs(
    explicit_dbs: Optional[Sequence[str]] = None,
    roots: Optional[Sequence[str]] = None,
) -> List[Tuple[str, str]]:
    """Resolve the set of structures.db files to operate on.

    Args:
        explicit_dbs: paths passed via --db (each must point to a structures.db).
        roots:        directory paths passed via --root; for each root we collect
                      every immediate subdirectory <root>/<formula>/structures.db,
                      and also <root>/structures.db if present.

    Returns:
        List of (formula_name, db_path) tuples. formula_name is the basename
        of the directory containing the DB (or the DB's own basename when the
        DB is given directly without a parent formula folder).
    """
    found: List[Tuple[str, str]] = []
    seen = set()

    def _add(formula: str, db: str) -> None:
        ap = os.path.abspath(db)
        if ap in seen:
            return
        seen.add(ap)
        found.append((formula, ap))

    for db in explicit_dbs or ():
        if not os.path.isfile(db):
            raise FileNotFoundError(f"--db path is not a file: {db}")
        formula = os.path.basename(os.path.dirname(os.path.abspath(db))) or "<db>"
        _add(formula, db)

    for root in roots or ():
        if not os.path.isdir(root):
            raise FileNotFoundError(f"--root is not a directory: {root}")
        direct = os.path.join(root, DB_NAME)
        if os.path.isfile(direct):
            _add(os.path.basename(os.path.abspath(root)) or "<root>", direct)
        for entry in sorted(os.listdir(root)):
            sub = os.path.join(root, entry)
            if not os.path.isdir(sub):
                continue
            db = os.path.join(sub, DB_NAME)
            if os.path.isfile(db):
                _add(entry, db)

    return found

def _parse_sg_spec(spec: str) -> List[Tuple[int, int]]:
    """Parse '225,227,100-150' into list of inclusive (lo, hi) ranges."""
    out: List[Tuple[int, int]] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            lo, hi = int(a), int(b)
            if lo > hi:
                lo, hi = hi, lo
            out.append((lo, hi))
        else:
            v = int(part)
            out.append((v, v))
    return out


@dataclass
class FilterSpec:
    """User-supplied filter for SELECT queries.

    All fields are optional; an empty spec selects every row.
    """
    formula_glob: Optional[str] = None        
    sg_spec: Optional[str] = None            
    name_glob: Optional[str] = None          
    stages: Optional[List[str]] = None        
    energy_min: Optional[float] = None
    energy_max: Optional[float] = None
    epa_min: Optional[float] = None
    epa_max: Optional[float] = None
    order_by: Optional[str] = None            
    limit: Optional[int] = None
    offset: int = 0

    def matches_formula(self, formula_name: str) -> bool:
        if not self.formula_glob:
            return True
        return fnmatch.fnmatch(formula_name, self.formula_glob)

    def to_sql(self, available_cols: Sequence[str]) -> Tuple[str, list]:
        """Return (where_clause, params). Where clause is empty string if no filter."""
        clauses: List[str] = []
        params: list = []

        if self.sg_spec:
            ranges = _parse_sg_spec(self.sg_spec)
            if ranges:
                ors = []
                for lo, hi in ranges:
                    if lo == hi:
                        ors.append("sg_number = ?")
                        params.append(lo)
                    else:
                        ors.append("sg_number BETWEEN ? AND ?")
                        params.extend([lo, hi])
                clauses.append("(" + " OR ".join(ors) + ")")

        if self.name_glob:
            clauses.append("cif_name GLOB ?")
            params.append(self.name_glob)

        if self.stages:
            placeholders = ",".join("?" * len(self.stages))
            clauses.append(f"stage IN ({placeholders})")
            params.extend(self.stages)

        if self.energy_min is not None:
            clauses.append("energy >= ?")
            params.append(self.energy_min)
        if self.energy_max is not None:
            clauses.append("energy <= ?")
            params.append(self.energy_max)
        if self.epa_min is not None:
            clauses.append("energy_per_atom >= ?")
            params.append(self.epa_min)
        if self.epa_max is not None:
            clauses.append("energy_per_atom <= ?")
            params.append(self.epa_max)

        where = ("WHERE " + " AND ".join(clauses)) if clauses else ""
        return where, params

    def order_clause(self, available_cols: Sequence[str]) -> str:
        if not self.order_by:
            return "ORDER BY sg_number, cif_name"
        col, _, direction = self.order_by.partition(":")
        col = _resolve_column(col.strip())
        direction = (direction or "asc").strip().lower()
        if col not in available_cols:
            raise ValueError(
                f"--order-by column '{col}' not present in DB "
                f"(available: {', '.join(available_cols)})"
            )
        if direction not in ("asc", "desc"):
            raise ValueError("order direction must be 'asc' or 'desc'")
        return f"ORDER BY {col} {direction.upper()}"

    def limit_clause(self) -> str:
        if self.limit is None:
            return ""
        return f"LIMIT {int(self.limit)} OFFSET {int(self.offset)}"

def _select_columns_for(conn: sqlite3.Connection, requested: Optional[Sequence[str]]) -> List[str]:
    """Resolve which columns to project.

    Filters out columns that don't exist on this particular DB so HP-only
    columns are silently skipped on non-HP databases.
    """
    available = set(list_columns(conn))
    if requested:
        cols = []
        for raw in requested:
            c = _resolve_column(raw)
            if c in available:
                cols.append(c)
        for ident in ("sg_number", "cif_name"):
            if ident not in cols:
                cols.insert(0, ident)
        return cols
    cols = [c for c in CORE_COLUMNS if c in available]
    for c in OPTIONAL_COLUMNS:
        if c in available:
            cols.append(c)
    return cols


def iter_rows(
    dbs: Sequence[Tuple[str, str]],
    spec: FilterSpec,
    columns: Optional[Sequence[str]] = None,
    include_cif_content: bool = False,
    cif_source: str = "relaxed",
) -> Iterator[dict]:
    """Yield filtered rows across all given DBs as plain dicts.

    Each yielded dict is augmented with:
      formula_name -- the parent directory name (a.k.a. formula folder)
      db_path      -- absolute path to the source DB

    cif_source controls which CIF text is returned in the 'cif_content' key
    when include_cif_content=True:
      'relaxed' (default) -- the live cif_content column
      'initial'           -- the cif_content_initial snapshot if present;
                              falls back to cif_content for rows that have
                              never been relaxed (their live column still
                              holds the initial CIF).
    """
    if cif_source not in ("relaxed", "initial"):
        raise ValueError("cif_source must be 'relaxed' or 'initial'")
    for formula_name, db_path in dbs:
        if not spec.matches_formula(formula_name):
            continue
        conn = open_ro(db_path)
        try:
            available = list_columns(conn)
            cols = _select_columns_for(conn, columns)
            if include_cif_content:
                if cif_source == "initial" and "cif_content_initial" in available:
                    extra = ("COALESCE(cif_content_initial, cif_content) AS cif_content",)
                else:
                    extra = ("cif_content",)
                cols = [c for c in cols if c not in ("cif_content",)] + list(extra)
            where, params = spec.to_sql(available)
            order = spec.order_clause(available)
            limit = spec.limit_clause()
            sql = f"SELECT {', '.join(cols)} FROM structures {where} {order} {limit}".strip()
            for row in conn.execute(sql, params):
                rec = {k: row[k] for k in row.keys()}
                rec["formula_name"] = formula_name
                rec["db_path"] = db_path
                yield rec
        finally:
            conn.close()


def count_rows(dbs: Sequence[Tuple[str, str]], spec: FilterSpec) -> int:
    """Return total matching rows across all DBs (ignores LIMIT/OFFSET)."""
    total = 0
    spec_no_limit = FilterSpec(
        formula_glob=spec.formula_glob,
        sg_spec=spec.sg_spec,
        name_glob=spec.name_glob,
        stages=spec.stages,
        energy_min=spec.energy_min,
        energy_max=spec.energy_max,
        epa_min=spec.epa_min,
        epa_max=spec.epa_max,
    )
    for formula_name, db_path in dbs:
        if not spec_no_limit.matches_formula(formula_name):
            continue
        conn = open_ro(db_path)
        try:
            available = list_columns(conn)
            where, params = spec_no_limit.to_sql(available)
            sql = f"SELECT COUNT(*) FROM structures {where}"
            total += conn.execute(sql, params).fetchone()[0]
        finally:
            conn.close()
    return total
