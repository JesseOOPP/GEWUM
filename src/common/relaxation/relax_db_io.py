"""
GEWUM DB-driven relaxation I/O helpers.

Two CLI entry points used by the shell drivers (relax_umlip.sh etc.):

  worker: read one (sg, cif_name) task directly from structures.db,
          relax it in memory, then dump the result to a per-task pkl in
          tmp_dir. No disk *_relaxed.cif is written; all downstream
          stages read CIF content from the DB.

  commit: aggregate all pkls in tmp_dir, parse the chemical formula
          from each relaxed CIF (pymatgen), batch-update structures.db
          (stage='relaxed' + cif_content + energy + formula + HP fields),
          append rows to the per-formula energy_results.csv where the
          Relaxed_CIF_Path column carries a 'db://<sg>/<cif_name>' URI
          instead of a filesystem path, then optionally clean up.

This split avoids SQLite write contention: the parallel workers only
read from the DB (safe under WAL), and a single-process commit step
performs all UPDATEs in one transaction.
"""

import os
import sys
import pickle
import sqlite3
import argparse
import tempfile
import logging

from gewum.src.common.relaxation.umlip_relax import (
    optimize_from_content,
    append_energy_csv,
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def _formula_from_cif_text(cif_text):
    """Parse chemical formula from a CIF string using pymatgen. Returns '' on failure."""
    try:
        from pymatgen.core import Structure
    except Exception:
        return ''
    fd, tmp_path = tempfile.mkstemp(suffix='.cif')
    os.close(fd)
    try:
        with open(tmp_path, 'w', encoding='utf-8') as fh:
            fh.write(cif_text)
        s = Structure.from_file(tmp_path)
        return s.composition.formula.replace(' ', '')
    except Exception:
        return ''
    finally:
        try:
            os.remove(tmp_path)
        except OSError:
            pass

def _read_cif_text(db_path, sg, cif_name):
    """Read the CIF content for one (sg, cif_name) pair. Read-only, WAL-safe."""
    conn = sqlite3.connect(db_path)
    try:
        row = conn.execute(
            "SELECT cif_content FROM structures WHERE sg_number=? AND cif_name=?",
            (sg, cif_name),
        ).fetchone()
    finally:
        conn.close()
    return row[0] if row else None


def run_worker(args):
    """Process exactly one relaxation task and dump the result to a pkl."""
    cif_text = _read_cif_text(args.db, args.sg, args.name)
    if cif_text is None:
        print(f"[relax_worker] not found in DB: sg={args.sg} name={args.name}",
              file=sys.stderr)
        return 2

    base_name = args.name[:-4] if args.name.lower().endswith('.cif') else args.name

    try:
        relaxed_text, energy, energy_per_atom, hp_data = optimize_from_content(
            cif_text, base_name,
            mode=args.mode, fmax=args.fmax,
            max_steps=args.max_steps, pressure=args.pressure,
        )
    except Exception as exc:  
        print(f"[relax_worker] FAILED sg={args.sg} name={args.name}: {exc}",
              file=sys.stderr)
        return 3

    os.makedirs(args.tmp_dir, exist_ok=True)
    safe_name = args.name.replace('/', '_').replace('\\', '_')
    pkl_path = os.path.join(args.tmp_dir, f"{args.sg}__{safe_name}.pkl")
    with open(pkl_path, 'wb') as fh:
        pickle.dump({
            'sg': args.sg,
            'name': args.name,
            'base_name': base_name,
            'relaxed_cif': relaxed_text,
            'energy': energy,
            'energy_per_atom': energy_per_atom,
            'hp_data': hp_data,
            'mode': args.mode,
        }, fh)

    logging.info(f"[worker] done sg={args.sg} {args.name} -> {pkl_path}")
    return 0

def run_commit(args):
    """Aggregate all pkls in tmp_dir, batch-update DB, append CSV."""
    if not os.path.isdir(args.tmp_dir):
        logging.info(f"[commit] no tmp_dir: {args.tmp_dir} (nothing to do)")
        return 0

    pkl_files = sorted(
        os.path.join(args.tmp_dir, f)
        for f in os.listdir(args.tmp_dir)
        if f.endswith('.pkl')
    )
    if not pkl_files:
        logging.info(f"[commit] no pkls in {args.tmp_dir}")
        return 0

    results = []
    for p in pkl_files:
        try:
            with open(p, 'rb') as fh:
                results.append(pickle.load(fh))
        except Exception as exc:  # noqa: BLE001
            logging.warning(f"[commit] skip unreadable pkl {p}: {exc}")

    if not results:
        return 0

    from gewum.src.common.cif_db import CifDatabase

    formula_dir = os.path.dirname(os.path.abspath(args.db))
    records = []
    for r in results:
        hp = r.get('hp_data') or {}
        records.append({
            'sg': r['sg'],
            'name': r['name'],
            'cif_content': r['relaxed_cif'],
            'energy': r['energy'],
            'energy_per_atom': r['energy_per_atom'],
            'formula': _formula_from_cif_text(r['relaxed_cif']),
            'final_pressure': hp.get('final_pressure'),
            'enthalpy_per_atom': hp.get('enthalpy_per_atom'),
            'corrected_enthalpy_per_atom': hp.get('corrected_enthalpy_per_atom'),
        })
    with CifDatabase(formula_dir) as db:
        db.update_relaxed_batch(records)
    logging.info(f"[commit] DB updated: {len(records)} records -> {args.db}")

    if args.csv:
        for r in results:
            db_uri = f"db://{r['sg']}/{r['name']}"
            append_energy_csv(
                args.csv, r['mode'], r['base_name'],
                r['energy'], r['energy_per_atom'],
                r['hp_data'], db_uri,
            )
        logging.info(f"[commit] CSV appended: {len(results)} rows -> {args.csv}")

    if args.cleanup:
        for p in pkl_files:
            try:
                os.remove(p)
            except OSError:
                pass
        try:
            os.rmdir(args.tmp_dir)
        except OSError:
            pass

    return 0

def main():
    parser = argparse.ArgumentParser(description='GEWUM DB-driven relaxation I/O')
    sub = parser.add_subparsers(dest='command', required=True)

    pw = sub.add_parser('worker', help='Run one relaxation task from DB')
    pw.add_argument('--db', required=True, help='Path to structures.db')
    pw.add_argument('--sg', type=int, required=True)
    pw.add_argument('--name', required=True, help='cif_name primary key value')
    pw.add_argument('--tmp-dir', required=True, help='Per-formula tmp dir for result pkls')
    pw.add_argument('--mode', type=int, default=2, choices=[1, 2, 3])
    pw.add_argument('--fmax', type=float, default=0.05)
    pw.add_argument('--max-steps', type=int, default=200)
    pw.add_argument('--pressure', type=float, default=0.0)

    pc = sub.add_parser('commit', help='Aggregate pkls and batch-update DB + CSV')
    pc.add_argument('--db', required=True, help='Path to structures.db')
    pc.add_argument('--tmp-dir', required=True)
    pc.add_argument('--csv', default='', help='Path to energy_results.csv (empty=skip CSV)')
    pc.add_argument('--cleanup', action='store_true', help='Remove pkls and tmp_dir on success')

    args = parser.parse_args()
    if args.command == 'worker':
        sys.exit(run_worker(args))
    elif args.command == 'commit':
        sys.exit(run_commit(args))
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
