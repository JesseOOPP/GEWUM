"""Batch CIF exporter for GEWUM DB admin (P1: read-only producer).

Streams cif_content from one or more structures.db files matching a
FilterSpec, then writes them either as plain files in a directory tree
or as a single zip archive. Optionally produces a manifest CSV indexing
every emitted file back to its source row.

Design notes:
  - Source DBs are opened read-only (see query.open_ro). This module
    never mutates user data.
  - Memory: rows are streamed via a cursor; cif_content is a TEXT column
    which is read on demand per row.
  - dry-run prints what would be written without touching the filesystem.
"""

from __future__ import annotations

import csv
import os
import re
import sys
import zipfile
from concurrent.futures import ThreadPoolExecutor
from typing import Iterable, List, Optional, Sequence, Tuple

from .query import FilterSpec, count_rows, iter_rows

LAYOUTS = ("flat", "by-formula", "by-sg", "by-stage", "nested")
OVERWRITE_MODES = ("skip", "replace", "fail")

_BAD_CHARS = re.compile(r'[\\/:*?"<>|\x00-\x1f]')


def _sanitize(name: str) -> str:
    """Strip path separators and control chars; collapse whitespace."""
    name = _BAD_CHARS.sub("_", str(name))
    name = name.strip().strip(".")
    return name or "_"


def _format_name(template: str, row: dict) -> str:
    """Apply str.format substitution with sanitised values.

    Recognised placeholders:
      {formula} {sg} {name} {stage} {epa} {energy}
    The original {name} keeps its '.cif' suffix if present; the template
    is responsible for re-adding '.cif' if it was stripped via custom format.
    """
    base_name = row.get("cif_name") or "structure.cif"
    epa = row.get("energy_per_atom")
    energy = row.get("energy")
    fields = {
        "formula": _sanitize(row.get("formula_name") or row.get("formula") or ""),
        "sg":      str(row.get("sg_number") or 0),
        "name":    _sanitize(os.path.splitext(base_name)[0]),
        "cif":     _sanitize(base_name), 
        "stage":   _sanitize(row.get("stage") or ""),
        "epa":     "NA" if epa is None else f"{epa:.6g}",
        "energy":  "NA" if energy is None else f"{energy:.6g}",
    }
    try:
        out = template.format(**fields)
    except KeyError as e:
        raise ValueError(f"Unknown placeholder {{{e.args[0]}}} in --name-template") from e
    if not out.lower().endswith(".cif"):
        out += ".cif"
    return _sanitize(out)


def _layout_subdir(layout: str, row: dict) -> str:
    if layout == "flat":
        return ""
    formula = _sanitize(row.get("formula_name") or row.get("formula") or "unknown")
    sg = str(row.get("sg_number") or 0)
    stage = _sanitize(row.get("stage") or "unknown")
    if layout == "by-formula":
        return formula
    if layout == "by-sg":
        return f"sg_{sg}"
    if layout == "by-stage":
        return stage
    if layout == "nested":
        return os.path.join(formula, f"sg_{sg}")
    raise ValueError(f"Unknown layout: {layout}")


def _resolve_target(out_dir: str, row: dict, layout: str, template: str) -> str:
    sub = _layout_subdir(layout, row)
    name = _format_name(template, row)
    return os.path.join(out_dir, sub, name) if sub else os.path.join(out_dir, name)

MANIFEST_COLUMNS = (
    "target_path", "formula_name", "sg_number", "cif_name",
    "stage", "energy_per_atom", "energy", "db_path",
)


def _open_manifest(path: Optional[str]):
    if not path:
        return None, None
    fp = open(path, "w", encoding="utf-8", newline="")
    writer = csv.writer(fp, lineterminator="\n")
    writer.writerow(MANIFEST_COLUMNS)
    return fp, writer


def _record_manifest(writer, target: str, row: dict) -> None:
    if writer is None:
        return
    writer.writerow([
        target,
        row.get("formula_name") or "",
        row.get("sg_number") or "",
        row.get("cif_name") or "",
        row.get("stage") or "",
        row.get("energy_per_atom") if row.get("energy_per_atom") is not None else "",
        row.get("energy") if row.get("energy") is not None else "",
        row.get("db_path") or "",
    ])

def export(
    dbs: Sequence[Tuple[str, str]],
    spec: FilterSpec,
    out_dir: Optional[str] = None,
    zip_path: Optional[str] = None,
    layout: str = "by-formula",
    name_template: str = "{name}.cif",
    workers: int = 1,
    overwrite: str = "skip",
    manifest: Optional[str] = None,
    max_files: int = 5000,
    confirm_above_max: bool = False,
    dry_run: bool = False,
    source: str = "relaxed",
    log=print,
) -> dict:
    """Export filtered CIFs to disk or a zip archive.

    Returns a small summary dict suitable for printing at the end.
    """
    if (out_dir is None) == (zip_path is None):
        raise ValueError("exactly one of --out / --zip must be supplied")
    if layout not in LAYOUTS:
        raise ValueError(f"--layout must be one of {LAYOUTS}")
    if overwrite not in OVERWRITE_MODES:
        raise ValueError(f"--overwrite must be one of {OVERWRITE_MODES}")
    if source not in ("relaxed", "initial"):
        raise ValueError("--source must be 'relaxed' or 'initial'")

    total = count_rows(dbs, spec)
    log(f"[export] {total} row(s) match the filter.")
    if total == 0:
        return {"total": 0, "written": 0, "skipped": 0, "failed": 0}

    if total > max_files and not confirm_above_max:
        raise RuntimeError(
            f"Filter matches {total} files which exceeds --max-files={max_files}. "
            f"Re-run with --confirm to proceed, or tighten the filter."
        )

    if dry_run:
        log("[export] --dry-run: showing first 20 target paths, no files will be written.")
        shown = 0
        for row in iter_rows(dbs, spec, include_cif_content=False):
            target = _resolve_target(out_dir or ".", row, layout, name_template)
            log(f"  {target}")
            shown += 1
            if shown >= 20:
                break
        log(f"[export] dry-run done. Would write {total} file(s) (source={source}).")
        return {"total": total, "written": 0, "skipped": 0, "failed": 0, "dry_run": True}

    manifest_fp, manifest_writer = _open_manifest(manifest)
    written = skipped = failed = 0

    try:
        if zip_path:
            written, skipped, failed = _export_zip(
                dbs, spec, zip_path, layout, name_template,
                overwrite, manifest_writer, log, source,
            )
        else:
            written, skipped, failed = _export_dir(
                dbs, spec, out_dir, layout, name_template,
                overwrite, workers, manifest_writer, log, source,
            )
    finally:
        if manifest_fp is not None:
            manifest_fp.close()

    log(
        f"[export] done. written={written}  skipped={skipped}  failed={failed}  "
        f"target={'zip:' + zip_path if zip_path else out_dir}"
    )
    if manifest:
        log(f"[export] manifest written to {manifest}")
    return {
        "total": total, "written": written,
        "skipped": skipped, "failed": failed,
    }

def _export_dir(dbs, spec, out_dir, layout, template, overwrite,
                workers, manifest_writer, log, source="relaxed"):
    os.makedirs(out_dir, exist_ok=True)

    jobs: List[Tuple[str, str, dict]] = []
    for row in iter_rows(dbs, spec, include_cif_content=True, cif_source=source):
        target = _resolve_target(out_dir, row, layout, template)
        jobs.append((target, row.get("cif_content") or "", row))

    written = skipped = failed = 0

    def write_one(item):
        nonlocal written, skipped, failed
        target, content, row = item
        try:
            if os.path.exists(target):
                if overwrite == "skip":
                    skipped += 1
                    return
                if overwrite == "fail":
                    raise FileExistsError(target)
            os.makedirs(os.path.dirname(target), exist_ok=True)
            with open(target, "w", encoding="utf-8") as f:
                f.write(content)
            _record_manifest(manifest_writer, target, row)
            written += 1
        except Exception as e:
            failed += 1
            log(f"  [!] failed to write {target}: {e}")

    if workers and workers > 1:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(write_one, jobs))
    else:
        for j in jobs:
            write_one(j)
    return written, skipped, failed

def _export_zip(dbs, spec, zip_path, layout, template, overwrite,
                manifest_writer, log, source="relaxed"):
    written = skipped = failed = 0
    seen_names = set()
    mode = "a" if (overwrite != "fail" and os.path.exists(zip_path)) else "w"
    with zipfile.ZipFile(zip_path, mode, compression=zipfile.ZIP_DEFLATED) as zf:
        existing = set(zf.namelist()) if mode == "a" else set()
        for row in iter_rows(dbs, spec, include_cif_content=True, cif_source=source):
            arcname = _resolve_target("", row, layout, template).lstrip(os.sep + "/")
            arcname = arcname.replace(os.sep, "/")
            try:
                if arcname in existing or arcname in seen_names:
                    if overwrite == "skip":
                        skipped += 1
                        continue
                    if overwrite == "fail":
                        raise FileExistsError(f"duplicate zip entry: {arcname}")
                zf.writestr(arcname, row.get("cif_content") or "")
                seen_names.add(arcname)
                _record_manifest(manifest_writer, arcname, row)
                written += 1
            except Exception as e:
                failed += 1
                log(f"  [!] failed zip-write {arcname}: {e}")
    return written, skipped, failed
