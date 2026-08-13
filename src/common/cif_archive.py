"""
GEWUM CIF Archive Module
Provides ZIP-based archival for CIF files to reduce filesystem inode pressure.

Each formula directory is packed into a single structures.zip containing
all spacegroup subdirectories and their CIF/CSV files.

Also supports the new SQLite-based structures.db (see cif_db.py).
DB-aware functions check for structures.db first, then fall back to ZIP
or loose filesystem files.

CLI Usage:
    python -m gewum.src.common.cif_archive pack   [--base-dir .]
    python -m gewum.src.common.cif_archive unpack  [--base-dir .]
"""
import os
import zipfile
import sqlite3
import argparse
import logging
from concurrent.futures import ThreadPoolExecutor

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

ARCHIVE_NAME = "structures.zip"
DB_NAME = "structures.db"

TOP_LEVEL_SKIP = {'final_cifs', '.energy_temp'}


def _is_formula_dir(dirpath):
    """
    Check if a directory looks like a formula directory (contains
    spacegroup subdirs, a structures.zip, or a structures.db).
    Skip known output directories and top-level special directories.
    """
    basename = os.path.basename(dirpath)
    if basename in TOP_LEVEL_SKIP or basename.startswith('0_'):
        return False
    if not os.path.isdir(dirpath):
        return False
    return True


def pack_directory(formula_dir):
    """
    Pack all files in a formula directory into structures.zip,
    then remove the original files and empty subdirectories.

    Args:
        formula_dir: Path to a formula directory (e.g., Na2Cl2/)

    Returns:
        Number of files packed, or 0 if nothing to pack
    """
    zip_path = os.path.join(formula_dir, ARCHIVE_NAME)

    files_to_pack = []
    for root, dirs, files in os.walk(formula_dir):
        dirs[:] = [d for d in dirs if not d.startswith('0_')]
        for fname in files:
            if fname == ARCHIVE_NAME:
                continue
            if root == formula_dir and fname.startswith('0_'):
                continue
            full_path = os.path.join(root, fname)
            rel_path = os.path.relpath(full_path, formula_dir)
            files_to_pack.append((full_path, rel_path))

    if not files_to_pack:
        return 0

    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for full_path, rel_path in files_to_pack:
            zf.write(full_path, rel_path)

    try:
        with zipfile.ZipFile(zip_path, 'r') as zf:
            if zf.testzip() is not None:
                logging.error(f"ZIP verification failed for {zip_path}, keeping originals")
                os.remove(zip_path)
                return 0
    except zipfile.BadZipFile:
        logging.error(f"Corrupted ZIP file {zip_path}, keeping originals")
        os.remove(zip_path)
        return 0

    def _remove_file(path):
        try:
            os.remove(path)
        except OSError:
            pass

    with ThreadPoolExecutor(max_workers=32) as executor:
        executor.map(_remove_file, [full_path for full_path, _ in files_to_pack])

    for root, dirs, files in os.walk(formula_dir, topdown=False):
        if root == formula_dir:
            continue
        try:
            if not os.listdir(root):
                os.rmdir(root)
        except OSError:
            pass

    logging.info(f"Packed {len(files_to_pack)} files into {zip_path}")
    return len(files_to_pack)


def unpack_directory(formula_dir):
    """
    Unpack structures.zip in a formula directory, then remove the zip.
    Uses parallel file writing to accelerate on shared filesystems (Lustre/GPFS).

    Args:
        formula_dir: Path to a formula directory

    Returns:
        Number of files unpacked, or 0 if no archive found
    """
    zip_path = os.path.join(formula_dir, ARCHIVE_NAME)

    if not os.path.isfile(zip_path):
        return 0

    file_data = {}
    with zipfile.ZipFile(zip_path, 'r') as zf:
        for member in zf.namelist():
            file_data[member] = zf.read(member)

    dirs_needed = set()
    for member in file_data:
        dest_path = os.path.join(formula_dir, member)
        parent_dir = os.path.dirname(dest_path)
        if parent_dir != formula_dir:
            dirs_needed.add(parent_dir)

    for d in sorted(dirs_needed):
        os.makedirs(d, exist_ok=True)

    def _write_file(member):
        dest_path = os.path.join(formula_dir, member)
        try:
            with open(dest_path, 'wb') as f:
                f.write(file_data[member])
        except OSError:
            pass

    with ThreadPoolExecutor(max_workers=32) as executor:
        executor.map(_write_file, file_data.keys())

    os.remove(zip_path)
    logging.info(f"Unpacked {len(file_data)} files from {zip_path}")
    return len(file_data)


def pack_all(base_dir='.'):
    """
    Pack all formula directories under base_dir.

    Args:
        base_dir: Working directory containing formula subdirectories
    """
    base_dir = os.path.abspath(base_dir)
    total_files = 0
    packed_dirs = 0

    for entry in sorted(os.listdir(base_dir)):
        dirpath = os.path.join(base_dir, entry)
        if not _is_formula_dir(dirpath):
            continue

        count = pack_directory(dirpath)
        if count > 0:
            total_files += count
            packed_dirs += 1

    logging.info(f"Pack complete: {total_files} files in {packed_dirs} directories")
    return total_files


def unpack_all(base_dir='.'):
    """
    Unpack all formula directories under base_dir that contain structures.zip.

    Args:
        base_dir: Working directory containing formula subdirectories
    """
    base_dir = os.path.abspath(base_dir)
    total_files = 0
    unpacked_dirs = 0

    for entry in sorted(os.listdir(base_dir)):
        dirpath = os.path.join(base_dir, entry)
        if not os.path.isdir(dirpath):
            continue

        zip_path = os.path.join(dirpath, ARCHIVE_NAME)
        if not os.path.isfile(zip_path):
            continue

        count = unpack_directory(dirpath)
        if count > 0:
            total_files += count
            unpacked_dirs += 1

    logging.info(f"Unpack complete: {total_files} files in {unpacked_dirs} directories")
    return total_files


def main():
    """CLI entry point for CIF archive operations."""
    parser = argparse.ArgumentParser(
        description='GEWUM CIF Archive - Pack/Unpack CIF files into ZIP archives',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'action',
        choices=['pack', 'unpack'],
        help='Action to perform:\n'
             '  pack   - Archive CIF files into ZIP per formula directory\n'
             '  unpack - Extract CIF files from ZIP archives'
    )
    parser.add_argument(
        '--base-dir', '-d',
        default='.',
        help='Base working directory (default: current directory)'
    )
    parser.add_argument(
        '--dir',
        default=None,
        help='Process a single directory instead of all formula dirs.\n'
             'When set, operates on this specific directory only.'
    )

    args = parser.parse_args()

    if args.dir:
        target = os.path.abspath(args.dir)
        if not os.path.isdir(target):
            logging.error(f"Directory not found: {target}")
            return
        if args.action == 'pack':
            pack_directory(target)
        elif args.action == 'unpack':
            unpack_directory(target)
    else:
        if args.action == 'pack':
            pack_all(args.base_dir)
        elif args.action == 'unpack':
            unpack_all(args.base_dir)

def list_cifs_in_zip(zip_path):
    """
    List all CIF member names inside a ZIP archive.

    Args:
        zip_path: Path to the ZIP file.

    Returns:
        Sorted list of member names ending with .cif
    """
    with zipfile.ZipFile(zip_path, 'r') as zf:
        return sorted(m for m in zf.namelist() if m.lower().endswith('.cif'))


def read_cif_from_zip(zip_path, member_name):
    """
    Read CIF content from a ZIP archive as a UTF-8 string.

    Args:
        zip_path: Path to the ZIP file.
        member_name: Member name inside the ZIP.

    Returns:
        CIF file content as string.
    """
    with zipfile.ZipFile(zip_path, 'r') as zf:
        return zf.read(member_name).decode('utf-8')


def entry_to_path(entry):
    """
    Convert a CIF entry to a virtual file path string.

    Args:
        entry: File path (str), 2-tuple (zip_path, member_name),
               3-tuple (db_path, sg_number, cif_name),
               or 4-tuple (db_path, sg_number, cif_name, column).
               When column == 'cif_content_initial' the path gets an
               '__initial__' marker dir to make it distinct from the
               relaxed counterpart of the same row.

    Returns:
        A path string usable for path-based logic (e.g., extracting space group from directory).
    """
    if isinstance(entry, tuple):
        if len(entry) == 4:
            db_path, sg_number, cif_name, column = entry
            base_dir = os.path.dirname(db_path)
            if column == 'cif_content_initial':
                return os.path.join(base_dir, '__initial__', str(sg_number), cif_name)
            return os.path.join(base_dir, str(sg_number), cif_name)
        if len(entry) == 3:
            db_path, sg_number, cif_name = entry
            return os.path.join(os.path.dirname(db_path), str(sg_number), cif_name)
        zip_path, member = entry
        return os.path.join(os.path.dirname(zip_path), member)
    return str(entry)


def entry_basename(entry):
    """
    Get the file basename from a CIF entry.

    Args:
        entry: File path (str), 2-tuple (zip_path, member_name),
               3-tuple (db_path, sg_number, cif_name),
               or 4-tuple (db_path, sg_number, cif_name, column).

    Returns:
        File basename string.
    """
    if isinstance(entry, tuple):
        if len(entry) == 4:
            return entry[2]
        if len(entry) == 3:
            return entry[2]
        return os.path.basename(entry[1])
    return os.path.basename(str(entry))


def load_structure(entry):
    """
    Load a pymatgen Structure from a CIF entry.

    Args:
        entry: File path (str), 2-tuple (zip_path, member_name),
               3-tuple (db_path, sg_number, cif_name),
               or 4-tuple (db_path, sg_number, cif_name, column)
               where column selects which DB column carries the CIF text
               (e.g. 'cif_content' or 'cif_content_initial').

    Returns:
        pymatgen Structure object.
    """
    from pymatgen.core.structure import Structure

    if isinstance(entry, tuple):
        if len(entry) == 4:
            db_path, sg_number, cif_name, column = entry
            content = read_cif_from_db(db_path, sg_number, cif_name, column=column)
            return Structure.from_str(content, fmt='cif')
        if len(entry) == 3:
            db_path, sg_number, cif_name = entry
            content = read_cif_from_db(db_path, sg_number, cif_name)
            return Structure.from_str(content, fmt='cif')
        zip_path, member = entry
        content = read_cif_from_zip(zip_path, member)
        return Structure.from_str(content, fmt='cif')
    else:
        return Structure.from_file(str(entry))


def find_cifs_zip_aware(directory, mode=None):
    """
    Find CIF files in a directory, transparently reading from structures.zip if present.

    Uses the same mode-based filtering as the original _find_cif_files functions.

    Args:
        directory: Root directory to search.
        mode: 'total', 'selected', 'relaxed', or None (all).

    Returns:
        Sorted list of entries: str paths for loose files, or (zip_path, member) tuples for ZIP.
    """
    from pathlib import Path as _Path

    EXCLUDE_DIRS = {'bond_mis'}
    zip_path = os.path.join(directory, ARCHIVE_NAME)

    if os.path.isfile(zip_path):
        entries = []
        with zipfile.ZipFile(zip_path, 'r') as zf:
            for member in zf.namelist():
                if not member.lower().endswith('.cif'):
                    continue
                parts = _Path(member).parts

                if EXCLUDE_DIRS & set(parts):
                    continue

                if mode == 'total':
                    if 'relaxed' in parts:
                        continue
                elif mode == 'selected':
                    if 'remove' in parts or 'relaxed' in parts:
                        continue
                elif mode == 'relaxed':
                    if 'relaxed' not in parts:
                        continue

                entries.append((zip_path, member))

        return sorted(entries, key=lambda e: e[1])
    else:
        cif_files = []
        for root, _, files in os.walk(directory):
            for f in files:
                if f.lower().endswith('.cif'):
                    full_path = os.path.join(root, f)
                    parts = _Path(full_path).parts

                    if EXCLUDE_DIRS & set(parts):
                        continue

                    if mode == 'total':
                        if 'relaxed' in parts:
                            continue
                    elif mode == 'selected':
                        if 'remove' in parts or 'relaxed' in parts:
                            continue
                    elif mode == 'relaxed':
                        if 'relaxed' not in parts:
                            continue

                    cif_files.append(full_path)

        return sorted(cif_files)


def detect_stages_zip_aware(cif_dir):
    """
    Detect workflow stages (remove/relaxed) in a directory, works with both
    loose files and structures.zip.

    Args:
        cif_dir: Path to a formula directory.

    Returns:
        dict with keys 'remove' and 'relaxed' (bool).
    """
    from pathlib import Path as _Path

    zip_path = os.path.join(cif_dir, ARCHIVE_NAME)
    has_remove = False
    has_relaxed = False

    if os.path.isfile(zip_path):
        with zipfile.ZipFile(zip_path, 'r') as zf:
            for member in zf.namelist():
                parts = _Path(member).parts
                if 'remove' in parts:
                    has_remove = True
                if 'relaxed' in parts:
                    has_relaxed = True
                if has_remove and has_relaxed:
                    break
    else:
        for root, dirs, _ in os.walk(cif_dir):
            if 'remove' in dirs:
                has_remove = True
            if 'relaxed' in dirs:
                has_relaxed = True
            if has_remove and has_relaxed:
                break

    return {'remove': has_remove, 'relaxed': has_relaxed}

def read_cif_from_db(db_path, sg_number, cif_name, column='cif_content'):
    """
    Read CIF content from a SQLite database as a UTF-8 string.

    Args:
        db_path: Path to the structures.db file.
        sg_number: Space group number.
        cif_name: CIF file name.
        column: Which CIF text column to read. Defaults to ``'cif_content'``
            (the live, possibly-relaxed structure). Use ``'cif_content_initial'``
            to read the pre-relax snapshot; falls back to ``cif_content`` via
            COALESCE for legacy rows that have no snapshot.

    Returns:
        CIF file content as string.
    """
    if column not in ('cif_content', 'cif_content_initial'):
        raise ValueError(
            f"read_cif_from_db: unsupported column {column!r}; "
            "expected 'cif_content' or 'cif_content_initial'"
        )
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        if column == 'cif_content_initial':
            cols = {r['name'] for r in conn.execute(
                "PRAGMA table_info(structures)"
            ).fetchall()}
            if 'cif_content_initial' in cols:
                sql = (
                    "SELECT COALESCE(cif_content_initial, cif_content) AS cif "
                    "FROM structures WHERE sg_number = ? AND cif_name = ?"
                )
            else:
                sql = (
                    "SELECT cif_content AS cif "
                    "FROM structures WHERE sg_number = ? AND cif_name = ?"
                )
        else:
            sql = (
                "SELECT cif_content AS cif "
                "FROM structures WHERE sg_number = ? AND cif_name = ?"
            )
        row = conn.execute(sql, (sg_number, cif_name)).fetchone()
        if row is None:
            raise KeyError(f"CIF not found in DB: SG={sg_number}, {cif_name}")
        return row['cif']
    finally:
        conn.close()


def find_cifs_db_aware(directory, mode=None):
    """
    Find CIF files in a directory, transparently reading from structures.db
    or structures.zip if present. Falls back to filesystem walk if neither exists.

    Uses the same mode-based filtering as _find_cif_files functions.

    Args:
        directory: Root directory to search.
        mode: 'total', 'selected', 'relaxed', or None (all).

    Returns:
        Sorted list of entries:
          - 3-tuple (db_path, sg_number, cif_name) for DB entries
          - 2-tuple (zip_path, member) for ZIP entries (backward compat)
          - str paths for loose files (backward compat)
    """
    db_path = os.path.join(directory, DB_NAME)
    if os.path.isfile(db_path):
        if mode == 'initial':
            conn = sqlite3.connect(db_path)
            conn.row_factory = sqlite3.Row
            try:
                cols = {r['name'] for r in conn.execute(
                    "PRAGMA table_info(structures)"
                ).fetchall()}
                if 'cif_content_initial' not in cols:
                    return []
                rows = conn.execute(
                    "SELECT sg_number, cif_name FROM structures "
                    "WHERE stage = 'relaxed' AND cif_content_initial IS NOT NULL "
                    "ORDER BY sg_number, cif_name"
                ).fetchall()
                return [(db_path, r['sg_number'], r['cif_name'], 'cif_content_initial')
                        for r in rows]
            finally:
                conn.close()

        if mode == 'selected':
            conn = sqlite3.connect(db_path)
            conn.row_factory = sqlite3.Row
            try:
                cols = {r['name'] for r in conn.execute(
                    "PRAGMA table_info(structures)"
                ).fetchall()}
                has_snapshot_col = 'cif_content_initial' in cols

                rows = conn.execute(
                    "SELECT sg_number, cif_name, stage "
                    "FROM structures WHERE stage != 'removed' "
                    "ORDER BY sg_number, cif_name"
                ).fetchall()

                entries = []
                for r in rows:
                    sg = r['sg_number']
                    name = r['cif_name']
                    if has_snapshot_col and r['stage'] == 'relaxed':
                        snap = conn.execute(
                            "SELECT 1 FROM structures "
                            "WHERE sg_number=? AND cif_name=? "
                            "AND cif_content_initial IS NOT NULL "
                            "LIMIT 1", (sg, name)
                        ).fetchone()
                        if snap:
                            entries.append((db_path, sg, name,
                                            'cif_content_initial'))
                        else:
                            entries.append((db_path, sg, name))
                    else:
                        entries.append((db_path, sg, name))
                return entries
            finally:
                conn.close()

        if mode == 'total':
            conn = sqlite3.connect(db_path)
            conn.row_factory = sqlite3.Row
            try:
                cols = {r['name'] for r in conn.execute(
                    "PRAGMA table_info(structures)"
                ).fetchall()}
                has_snapshot_col = 'cif_content_initial' in cols

                rows = conn.execute(
                    "SELECT sg_number, cif_name, stage "
                    "FROM structures ORDER BY sg_number, cif_name"
                ).fetchall()

                entries = []
                for r in rows:
                    sg = r['sg_number']
                    name = r['cif_name']
                    if has_snapshot_col and r['stage'] == 'relaxed':
                        snap = conn.execute(
                            "SELECT 1 FROM structures "
                            "WHERE sg_number=? AND cif_name=? "
                            "AND cif_content_initial IS NOT NULL "
                            "LIMIT 1", (sg, name)
                        ).fetchone()
                        if snap:
                            entries.append((db_path, sg, name,
                                            'cif_content_initial'))
                        else:
                            entries.append((db_path, sg, name))
                    else:
                        entries.append((db_path, sg, name))
                return entries
            finally:
                conn.close()

        stage_filter = {
            None:       None,
            'all':      None,
            'relaxed':  "stage = 'relaxed'",
        }.get(mode)

        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row
        try:
            if stage_filter:
                rows = conn.execute(
                    f"SELECT sg_number, cif_name FROM structures WHERE {stage_filter} ORDER BY sg_number, cif_name"
                ).fetchall()
            else:
                rows = conn.execute(
                    "SELECT sg_number, cif_name FROM structures ORDER BY sg_number, cif_name"
                ).fetchall()
            return [(db_path, r['sg_number'], r['cif_name']) for r in rows]
        finally:
            conn.close()

    return find_cifs_zip_aware(directory, mode)


def detect_stages_db_aware(cif_dir):
    """
    Detect workflow stages (remove/relaxed/initial) in a directory, works with
    structures.db, structures.zip, or loose files.

    Args:
        cif_dir: Path to a formula directory.

    Returns:
        dict with keys 'remove', 'relaxed' and 'initial' (bool).
        'initial' is True only when the DB has a non-empty cif_content_initial
        snapshot for at least one relaxed row (i.e. true pre-vs-post comparison
        data is available).
    """
    db_path = os.path.join(cif_dir, DB_NAME)
    if os.path.isfile(db_path):
        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row
        try:
            stages = {r['stage'] for r in conn.execute(
                "SELECT DISTINCT stage FROM structures"
            ).fetchall()}
            cols = {r['name'] for r in conn.execute(
                "PRAGMA table_info(structures)"
            ).fetchall()}
            has_initial_snapshot = False
            if 'cif_content_initial' in cols and 'relaxed' in stages:
                row = conn.execute(
                    "SELECT 1 FROM structures "
                    "WHERE stage = 'relaxed' AND cif_content_initial IS NOT NULL "
                    "LIMIT 1"
                ).fetchone()
                has_initial_snapshot = row is not None
            return {
                'remove':  'removed' in stages,
                'relaxed': 'relaxed' in stages,
                'initial': has_initial_snapshot,
            }
        finally:
            conn.close()

    info = detect_stages_zip_aware(cif_dir)
    info.setdefault('initial', False)
    return info

def migrate_zip_to_db(formula_dir, remove_zip=False):
    """
    Migrate a formula directory from structures.zip to structures.db.

    Args:
        formula_dir: Path to formula directory.
        remove_zip: If True, delete the ZIP after migration.

    Returns:
        Number of CIF files migrated.
    """
    from gewum.src.common.cif_db import CifDatabase
    db = CifDatabase(formula_dir)
    return db.migrate_from_zip(remove_zip=remove_zip)


if __name__ == '__main__':
    main()
