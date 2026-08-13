"""
GEWUM CIF Database Module
SQLite-based storage for CIF files, replacing ZIP-based structures.zip.

Each formula directory stores a single structures.db containing all
CIF files with stage tracking (initial / removed / relaxed).

Key advantages over ZIP:
  - No pack/unpack needed; all operations are in-place SQL queries
  - Atomic updates with transactions
  - Efficient filtering by stage and space group
  - No filesystem inode pressure while maintaining random access

CLI Usage:
    python -m gewum.src.common.cif_db init     --dir <formula_dir>
    python -m gewum.src.common.cif_db insert   --dir <formula_dir> --sg <N> --cif <file.cif>
    python -m gewum.src.common.cif_db query    --dir <formula_dir> [--stage <stage>] [--sg <N>]
    python -m gewum.src.common.cif_db update   --dir <formula_dir> --sg <N> --cif <name> --stage <stage>
    python -m gewum.src.common.cif_db migrate  --dir <formula_dir>  # from structures.zip
    python -m gewum.src.common.cif_db migrate-all --base-dir <dir>   # batch migration
    python -m gewum.src.common.cif_db stats    --dir <formula_dir>
"""

import os
import sys
import sqlite3
import argparse
import logging
import zipfile
from concurrent.futures import ThreadPoolExecutor

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

DB_NAME = "structures.db"
ARCHIVE_NAME = "structures.zip" 

TOP_LEVEL_SKIP = {'final_cifs', '.energy_temp'}


class CifDatabase:
    """SQLite-backed CIF storage for a single formula directory."""

    def __init__(self, formula_dir):
        """
        Args:
            formula_dir: Path to a formula directory (e.g., Na2Cl2/)
        """
        self.formula_dir = os.path.abspath(formula_dir)
        self.db_path = os.path.join(self.formula_dir, DB_NAME)
        self._conn = None

    @property
    def conn(self):
        if self._conn is None:
            self._conn = sqlite3.connect(self.db_path)
            self._conn.row_factory = sqlite3.Row
            self._conn.execute("PRAGMA journal_mode=WAL")
            self._conn.execute("PRAGMA synchronous=NORMAL")
        return self._conn

    def close(self):
        if self._conn:
            self._conn.close()
            self._conn = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    _OPTIONAL_COLUMNS = {
        'final_pressure':              'REAL',
        'enthalpy_per_atom':           'REAL',
        'corrected_enthalpy_per_atom': 'REAL',
        'cif_content_initial':         'TEXT',
    }

    def _ensure_columns(self, columns=None):
        """Add any missing columns to structures table (SQLite-safe idempotent ALTER)."""
        columns = columns or self._OPTIONAL_COLUMNS
        existing = {r['name'] for r in self.conn.execute(
            "PRAGMA table_info(structures)"
        ).fetchall()}
        for name, type_ in columns.items():
            if name not in existing:
                self.conn.execute(f"ALTER TABLE structures ADD COLUMN {name} {type_}")
        self.conn.commit()

    def init_db(self):
        """Create tables and indexes if they don't exist."""
        self.conn.executescript("""
            CREATE TABLE IF NOT EXISTS structures (
                sg_number   INTEGER NOT NULL,
                cif_name    TEXT    NOT NULL,
                cif_content TEXT    NOT NULL,
                stage       TEXT    NOT NULL DEFAULT 'initial',
                energy      REAL,
                energy_per_atom REAL,
                formula     TEXT,
                created_at  TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                PRIMARY KEY (sg_number, cif_name)
            );
            CREATE INDEX IF NOT EXISTS idx_structures_stage
                ON structures(stage);
            CREATE INDEX IF NOT EXISTS idx_structures_sg
                ON structures(sg_number);
            CREATE INDEX IF NOT EXISTS idx_structures_formula
                ON structures(formula);
        """)
        self.conn.commit()
        self._ensure_columns()

    def drop_db(self):
        """Delete the database file."""
        self.close()
        if os.path.isfile(self.db_path):
            os.remove(self.db_path)

    def insert(self, sg_number, cif_name, cif_content, stage='initial',
               energy=None, energy_per_atom=None, formula=None):
        """Insert or replace a CIF record."""
        self.conn.execute("""
            INSERT OR REPLACE INTO structures
                (sg_number, cif_name, cif_content, stage, energy, energy_per_atom, formula)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """, (sg_number, cif_name, cif_content, stage, energy, energy_per_atom, formula))
        self.conn.commit()

    def insert_batch(self, records):
        """Insert multiple records in a single transaction.
        
        Args:
            records: List of (sg_number, cif_name, cif_content, stage, energy, energy_per_atom, formula)
        """
        with self.conn:
            self.conn.executemany("""
                INSERT OR REPLACE INTO structures
                    (sg_number, cif_name, cif_content, stage, energy, energy_per_atom, formula)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, records)

    def get(self, sg_number, cif_name):
        """Get a single CIF record by primary key."""
        row = self.conn.execute(
            "SELECT * FROM structures WHERE sg_number = ? AND cif_name = ?",
            (sg_number, cif_name)
        ).fetchone()
        return dict(row) if row else None

    def get_content(self, sg_number, cif_name):
        """Get just the CIF content string."""
        row = self.conn.execute(
            "SELECT cif_content FROM structures WHERE sg_number = ? AND cif_name = ?",
            (sg_number, cif_name)
        ).fetchone()
        return row['cif_content'] if row else None

    def query(self, stage=None, sg_number=None):
        """Query CIF records with optional filters.
        
        Returns:
            List of dict rows.
        """
        conditions = []
        params = []
        if stage is not None:
            conditions.append("stage = ?")
            params.append(stage)
        if sg_number is not None:
            conditions.append("sg_number = ?")
            params.append(sg_number)

        if conditions:
            where = "WHERE " + " AND ".join(conditions)
        else:
            where = ""

        rows = self.conn.execute(
            f"SELECT * FROM structures {where} ORDER BY sg_number, cif_name",
            params
        ).fetchall()
        return [dict(r) for r in rows]

    def list_entries(self, mode=None):
        """List CIF entries with optional stage filtering.

        DB-era semantics (stage is a single label per row, not a directory copy):
          - None / 'all': every row, regardless of stage
          - 'total':      every row (DB has no 'redundant' relaxed copy to exclude)
          - 'selected':   rows that survived the selection step (stage != 'removed').
                          Includes 'initial' (waiting for relax), 'relaxed' (done) and
                          'bond_mis' (relax done but flagged). |selected| >= |relaxed|.
          - 'relaxed':    only stage='relaxed' (excludes 'bond_mis' and 'initial')

        Returns:
            List of (db_path, sg_number, cif_name) 3-tuples.
        """
        stage_filter = {
            None:       None,
            'all':      None,
            'total':    None,
            'selected': "stage != 'removed'",
            'relaxed':  "stage = 'relaxed'",
        }.get(mode)

        if stage_filter:
            rows = self.conn.execute(
                f"SELECT sg_number, cif_name FROM structures WHERE {stage_filter} ORDER BY sg_number, cif_name"
            ).fetchall()
        else:
            rows = self.conn.execute(
                "SELECT sg_number, cif_name FROM structures ORDER BY sg_number, cif_name"
            ).fetchall()

        return [(self.db_path, r['sg_number'], r['cif_name']) for r in rows]

    def count(self, stage=None):
        """Count CIF records, optionally filtered by stage."""
        if stage:
            row = self.conn.execute(
                "SELECT COUNT(*) as cnt FROM structures WHERE stage = ?", (stage,)
            ).fetchone()
        else:
            row = self.conn.execute(
                "SELECT COUNT(*) as cnt FROM structures"
            ).fetchone()
        return row['cnt']

    def update_stage(self, sg_number, cif_name, new_stage):
        """Change the stage of a CIF record."""
        self.conn.execute(
            "UPDATE structures SET stage = ? WHERE sg_number = ? AND cif_name = ?",
            (new_stage, sg_number, cif_name)
        )
        self.conn.commit()

    def update_stage_batch(self, items, new_stage):
        """Batch update stage for multiple (sg_number, cif_name) pairs."""
        with self.conn:
            self.conn.executemany(
                "UPDATE structures SET stage = ? WHERE sg_number = ? AND cif_name = ?",
                [(new_stage, sg, name) for sg, name in items]
            )

    def update_energy(self, sg_number, cif_name, energy, energy_per_atom=None):
        """Update energy data for a CIF record."""
        self.conn.execute(
            "UPDATE structures SET energy = ?, energy_per_atom = ? WHERE sg_number = ? AND cif_name = ?",
            (energy, energy_per_atom, sg_number, cif_name)
        )
        self.conn.commit()

    def update_relaxed(self, sg_number, cif_name, cif_content,
                       energy=None, energy_per_atom=None):
        """Write back relaxation result atomically and mark stage='relaxed'.

        Before overwriting cif_content with the relaxed structure, the
        existing (initial) cif_content is snapshotted into cif_content_initial
        via COALESCE so the pre-relax CIF text is preserved. COALESCE makes
        the snapshot a no-op on subsequent re-relaxations.
        """
        self._ensure_columns()
        self.conn.execute(
            "UPDATE structures SET "
            "cif_content_initial=COALESCE(cif_content_initial, cif_content), "
            "cif_content=?, energy=?, energy_per_atom=?, stage='relaxed' "
            "WHERE sg_number=? AND cif_name=?",
            (cif_content, energy, energy_per_atom, sg_number, cif_name)
        )
        self.conn.commit()

    def update_relaxed_batch(self, records):
        """Batch write back relaxation results in a single transaction.

        Args:
            records: iterable of dicts with required keys
                sg, name, cif_content, energy, energy_per_atom
            and optional keys
                formula, final_pressure, enthalpy_per_atom, corrected_enthalpy_per_atom

        Behaviour:
            - stage is always set to 'relaxed'
            - formula is overwritten only when provided (COALESCE keeps the
              existing DB value if the caller passed None)
            - HP columns are written verbatim (NULL when not provided)
        """
        self._ensure_columns()

        def _f(x):
            if x is None:
                return None
            if isinstance(x, (bytes, bytearray)):
                import struct
                try:
                    if len(x) == 4:
                        return float(struct.unpack('<f', x)[0])
                    if len(x) == 8:
                        return float(struct.unpack('<d', x)[0])
                except struct.error:
                    return None
                return None
            try:
                return float(x)
            except (TypeError, ValueError):
                return None

        rows = [
            (
                r['cif_content'], _f(r['energy']), _f(r['energy_per_atom']),
                r.get('formula'),
                _f(r.get('final_pressure')),
                _f(r.get('enthalpy_per_atom')),
                _f(r.get('corrected_enthalpy_per_atom')),
                r['sg'], r['name'],
            )
            for r in records
        ]
        with self.conn:
            self.conn.executemany(
                "UPDATE structures SET "
                "cif_content_initial=COALESCE(cif_content_initial, cif_content), "
                "cif_content=?, energy=?, energy_per_atom=?, "
                "formula=COALESCE(?, formula), "
                "final_pressure=?, enthalpy_per_atom=?, corrected_enthalpy_per_atom=?, "
                "stage='relaxed' "
                "WHERE sg_number=? AND cif_name=?",
                rows,
            )

    def query_initial_tasks(self, sg_number=None):
        """Return pending relaxation tasks (stage='initial') with their CIF content.

        Returns:
            List of (sg_number, cif_name, cif_content) tuples.
        """
        if sg_number is not None:
            rows = self.conn.execute(
                "SELECT sg_number, cif_name, cif_content FROM structures "
                "WHERE stage='initial' AND sg_number=? ORDER BY sg_number, cif_name",
                (sg_number,)
            ).fetchall()
        else:
            rows = self.conn.execute(
                "SELECT sg_number, cif_name, cif_content FROM structures "
                "WHERE stage='initial' ORDER BY sg_number, cif_name"
            ).fetchall()
        return [(r['sg_number'], r['cif_name'], r['cif_content']) for r in rows]

    def query_relaxed_for_post(self):
        """Return relaxed records (stage='relaxed') for energy post-processing.

        Returns:
            List of dicts with keys sg_number/cif_name/cif_content/energy/
            energy_per_atom/formula/final_pressure/enthalpy_per_atom/
            corrected_enthalpy_per_atom.
            Only entries with non-NULL energy_per_atom are returned.
        """
        self._ensure_columns()
        rows = self.conn.execute(
            "SELECT sg_number, cif_name, cif_content, energy, energy_per_atom, formula, "
            "final_pressure, enthalpy_per_atom, corrected_enthalpy_per_atom "
            "FROM structures WHERE stage='relaxed' AND energy_per_atom IS NOT NULL "
            "ORDER BY sg_number, cif_name"
        ).fetchall()
        return [dict(r) for r in rows]

    def delete(self, sg_number, cif_name):
        """Delete a single CIF record."""
        self.conn.execute(
            "DELETE FROM structures WHERE sg_number = ? AND cif_name = ?",
            (sg_number, cif_name)
        )
        self.conn.commit()

    def get_stages(self):
        """Get set of distinct stages in the database."""
        rows = self.conn.execute("SELECT DISTINCT stage FROM structures").fetchall()
        return {r['stage'] for r in rows}

    def stats(self):
        """Return summary statistics."""
        total = self.count()
        stages = {}
        for row in self.conn.execute(
            "SELECT stage, COUNT(*) as cnt FROM structures GROUP BY stage"
        ).fetchall():
            stages[row['stage']] = row['cnt']
        return {'total': total, 'stages': stages}

    def migrate_from_zip(self, remove_zip=False):
        """Import all CIF files from the existing structures.zip.
        
        Path convention inside ZIP:
            SG_number/file.cif          -> stage='initial'
            SG_number/relaxed/file.cif  -> stage='relaxed'
            SG_number/remove/file.cif   -> stage='removed'
        
        Args:
            remove_zip: If True, delete the ZIP after successful migration.
        
        Returns:
            Number of CIF files migrated.
        """
        zip_path = os.path.join(self.formula_dir, ARCHIVE_NAME)
        if not os.path.isfile(zip_path):
            logging.warning(f"No {ARCHIVE_NAME} found in {self.formula_dir}")
            return 0

        self.init_db()

        count = 0
        with zipfile.ZipFile(zip_path, 'r') as zf:
            for member in zf.namelist():
                if not member.lower().endswith('.cif'):
                    continue

                parts = member.replace('\\', '/').split('/')
                try:
                    sg_number = int(parts[0])
                except (ValueError, IndexError):
                    continue

                cif_name = parts[-1]
                stage = 'initial'
                if len(parts) >= 3 and parts[1] == 'relaxed':
                    stage = 'relaxed'
                elif len(parts) >= 3 and parts[1] == 'remove':
                    stage = 'removed'

                content = zf.read(member).decode('utf-8', errors='replace')
                self.insert(sg_number, cif_name, content, stage=stage)
                count += 1

        if count > 0 and remove_zip:
            os.remove(zip_path)
            logging.info(f"Migrated {count} CIFs and removed {zip_path}")
        else:
            logging.info(f"Migrated {count} CIFs from {zip_path}")

        return count

    def import_from_filesystem(self):
        """Import CIF files from the formula directory's subdirectories.
        
        Scans SG_number/*.cif and SG_number/relaxed/*.cif etc.
        Useful as an alternative to migrate_from_zip when files are loose.
        """
        self.init_db()
        count = 0

        for root, dirs, files in os.walk(self.formula_dir):
            dirs[:] = [d for d in dirs if not d.startswith('0_')]

            for fname in files:
                if not fname.lower().endswith('.cif'):
                    continue

                full_path = os.path.join(root, fname)
                rel_path = os.path.relpath(full_path, self.formula_dir)
                parts = rel_path.replace('\\', '/').split('/')

                try:
                    sg_number = int(parts[0])
                except (ValueError, IndexError):
                    continue

                stage = 'initial'
                if 'relaxed' in parts:
                    stage = 'relaxed'
                elif 'remove' in parts:
                    stage = 'removed'

                with open(full_path, 'r', encoding='utf-8', errors='replace') as f:
                    content = f.read()

                self.insert(sg_number, fname, content, stage=stage)
                count += 1

        logging.info(f"Imported {count} CIF files from filesystem into {self.db_path}")
        return count

def _is_formula_dir(dirpath):
    """Check if a directory looks like a formula directory."""
    basename = os.path.basename(dirpath)
    if basename in TOP_LEVEL_SKIP or basename.startswith('0_'):
        return False
    if not os.path.isdir(dirpath):
        return False
    return True


def find_formula_dirs(base_dir):
    """Find all formula directories under base_dir."""
    base_dir = os.path.abspath(base_dir)
    dirs = []
    for entry in sorted(os.listdir(base_dir)):
        dirpath = os.path.join(base_dir, entry)
        if _is_formula_dir(dirpath):
            dirs.append(dirpath)
    return dirs

def main():
    parser = argparse.ArgumentParser(
        description='GEWUM CIF Database - SQLite storage for CIF files',
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    p_init = subparsers.add_parser('init', help='Initialize structures.db in a formula directory')
    p_init.add_argument('--dir', required=True, help='Formula directory path')

    p_insert = subparsers.add_parser('insert', help='Insert a CIF file into the database')
    p_insert.add_argument('--dir', required=True, help='Formula directory path')
    p_insert.add_argument('--sg', type=int, required=True, help='Space group number')
    p_insert.add_argument('--cif', required=True, help='Path to CIF file to insert')
    p_insert.add_argument('--stage', default='initial', help='Stage (default: initial)')

    p_query = subparsers.add_parser('query', help='Query CIF records')
    p_query.add_argument('--dir', required=True, help='Formula directory path')
    p_query.add_argument('--stage', default=None, help='Filter by stage')
    p_query.add_argument('--sg', type=int, default=None, help='Filter by space group')

    p_update = subparsers.add_parser('update', help='Update a CIF record')
    p_update.add_argument('--dir', required=True, help='Formula directory path')
    p_update.add_argument('--sg', type=int, required=True, help='Space group number')
    p_update.add_argument('--cif', required=True, help='CIF file name')
    p_update.add_argument('--stage', default=None, help='New stage value')

    p_migrate = subparsers.add_parser('migrate', help='Migrate single formula dir from structures.zip to structures.db')
    p_migrate.add_argument('--dir', required=True, help='Formula directory path')
    p_migrate.add_argument('--remove-zip', action='store_true', help='Remove ZIP after migration')

    p_migrate_all = subparsers.add_parser('migrate-all', help='Batch migrate all formula dirs')
    p_migrate_all.add_argument('--base-dir', default='.', help='Base working directory')
    p_migrate_all.add_argument('--remove-zip', action='store_true', help='Remove ZIP after migration')

    p_stats = subparsers.add_parser('stats', help='Show database statistics')
    p_stats.add_argument('--dir', required=True, help='Formula directory path')

    p_import = subparsers.add_parser('import-fs', help='Import CIFs from filesystem into DB')
    p_import.add_argument('--dir', required=True, help='Formula directory path')

    args = parser.parse_args()

    if args.command == 'init':
        db = CifDatabase(args.dir)
        db.init_db()
        print(f"Initialized {db.db_path}")

    elif args.command == 'insert':
        db = CifDatabase(args.dir)
        db.init_db()
        cif_path = args.cif
        cif_name = os.path.basename(cif_path)
        with open(cif_path, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()
        db.insert(args.sg, cif_name, content, stage=args.stage)
        print(f"Inserted {cif_name} (SG={args.sg}, stage={args.stage})")

    elif args.command == 'query':
        db = CifDatabase(args.dir)
        rows = db.query(stage=args.stage, sg_number=args.sg)
        for r in rows:
            print(f"  SG={r['sg_number']:4d}  {r['cif_name']:40s}  stage={r['stage']}")

    elif args.command == 'update':
        db = CifDatabase(args.dir)
        db.update_stage(args.sg, args.cif, args.stage)
        print(f"Updated {args.cif} (SG={args.sg}) stage -> {args.stage}")

    elif args.command == 'migrate':
        db = CifDatabase(args.dir)
        count = db.migrate_from_zip(remove_zip=args.remove_zip)
        print(f"Migrated {count} CIFs from {args.dir}")

    elif args.command == 'migrate-all':
        base_dir = os.path.abspath(args.base_dir)
        formula_dirs = find_formula_dirs(base_dir)
        total = 0
        for formula_dir in formula_dirs:
            db = CifDatabase(formula_dir)
            count = db.migrate_from_zip(remove_zip=args.remove_zip)
            total += count
        print(f"Total: {total} CIFs migrated from {len(formula_dirs)} formula directories")

    elif args.command == 'stats':
        db = CifDatabase(args.dir)
        s = db.stats()
        print(f"Database: {db.db_path}")
        print(f"  Total CIFs: {s['total']}")
        for stage, cnt in sorted(s['stages'].items()):
            print(f"    {stage}: {cnt}")

    elif args.command == 'import-fs':
        db = CifDatabase(args.dir)
        count = db.import_from_filesystem()
        print(f"Imported {count} CIFs")

    else:
        parser.print_help()


if __name__ == '__main__':
    main()
