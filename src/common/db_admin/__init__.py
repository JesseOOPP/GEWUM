"""GEWUM DB Admin (P1: read-only + export).

User-facing tooling to inspect and selectively export CIFs from the
SQLite-backed structures.db files produced by RD/COMP workflows.

This package is intentionally read-only:
  - All sqlite connections are opened with mode=ro (URI form), so it is
    physically impossible for any sub-command here to mutate the DB.
  - There are no INSERT / UPDATE / DELETE / ALTER statements in this package.

Sub-modules:
  query     -- read-only connection factory, multi-DB discovery, FilterSpec -> SQL
  formatter -- table / csv / json / tsv renderers for query results
  stats     -- one/two-dimensional GROUP BY aggregation (pivot table)
  exporter  -- batch CIF export with layout/template/zip/parallel/manifest
"""

from .query import (
    FilterSpec,
    discover_dbs,
    open_ro,
    iter_rows,
    count_rows,
)
from .formatter import render
from .stats import aggregate
from .exporter import export

__all__ = [
    "FilterSpec",
    "discover_dbs",
    "open_ro",
    "iter_rows",
    "count_rows",
    "render",
    "aggregate",
    "export",
]
