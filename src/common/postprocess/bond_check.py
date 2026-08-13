"""
GEWUM Bond Length Check Module (DB-driven, stage 2).

Scans every structures.db under the given base directory, parses the
relaxed CIFs in memory, and marks structures whose minimum bond length
is below the threshold by setting their stage to 'bond_mis' (instead of
the legacy disk shutil.move into a bond_mis/ subfolder).

Each formula directory has its own structures.db, so per-formula work
can run in parallel processes without SQLite write contention.
"""

import os
import argparse
import warnings
import tempfile
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

from pymatgen.io.cif import CifParser

from gewum.src.common.cif_db import CifDatabase

warnings.filterwarnings(
    "ignore",
    message="Issues encountered while parsing CIF: .* fractional coordinates rounded to ideal values to avoid issues with finite precision.",
    category=UserWarning,
)


def calculate_min_bond_length(structure, search_radius=3.5):
    """Return the minimum pairwise distance in a pymatgen Structure."""
    min_distance = float('inf')
    for site in structure:
        for neighbor in structure.get_neighbors(site, r=search_radius):
            d = site.distance(neighbor[0])
            if d < min_distance:
                min_distance = d
    return min_distance


def _parse_cif_text(cif_text):
    """Parse a CIF string with pymatgen via a temp file."""
    fd, tmp_path = tempfile.mkstemp(suffix='.cif')
    os.close(fd)
    try:
        with open(tmp_path, 'w', encoding='utf-8') as fh:
            fh.write(cif_text)
        parser = CifParser(tmp_path)
        structures = parser.parse_structures()
        return structures[0] if structures else None
    except Exception:
        return None
    finally:
        try:
            os.remove(tmp_path)
        except OSError:
            pass


def process_formula_db(formula_dir, threshold):
    """Mark all stage='relaxed' entries in formula_dir/structures.db whose
    minimum bond length is below `threshold` as stage='bond_mis'.

    Returns:
        (n_checked, n_marked)
    """
    db_path = os.path.join(formula_dir, 'structures.db')
    if not os.path.isfile(db_path):
        return 0, 0

    print(f"[bond_check] {db_path}")
    bond_mis_items = []
    n_checked = 0

    with CifDatabase(formula_dir) as db:
        rows = db.query(stage='relaxed')
        for r in rows:
            n_checked += 1
            structure = _parse_cif_text(r['cif_content'])
            if structure is None:
                print(f"  parse error: SG={r['sg_number']} {r['cif_name']}")
                continue
            try:
                min_d = calculate_min_bond_length(structure)
            except Exception as exc:  # noqa: BLE001
                print(f"  neighbor error: SG={r['sg_number']} {r['cif_name']}: {exc}")
                continue
            if min_d < threshold:
                bond_mis_items.append((r['sg_number'], r['cif_name']))
                print(f"  bond_mis: SG={r['sg_number']} {r['cif_name']} "
                      f"min_bond={min_d:.3f}")

        if bond_mis_items:
            db.update_stage_batch(bond_mis_items, 'bond_mis')

    return n_checked, len(bond_mis_items)


def find_formula_dirs_with_db(base_dir):
    """Find every directory containing structures.db under base_dir.

    Accepts both layouts:
      base_dir/<formula>/structures.db   (multi-formula run)
      base_dir/structures.db             (single formula)
    """
    base_dir = os.path.abspath(base_dir)
    result = []
    if os.path.isfile(os.path.join(base_dir, 'structures.db')):
        result.append(base_dir)
    if os.path.isdir(base_dir):
        for entry in sorted(os.listdir(base_dir)):
            p = os.path.join(base_dir, entry)
            if os.path.isdir(p) and os.path.isfile(os.path.join(p, 'structures.db')):
                if p not in result:
                    result.append(p)
    return result


def main():
    parser = argparse.ArgumentParser(description='GEWUM Bond Length Check (DB-driven)')
    parser.add_argument('--threshold', '-t', type=float, default=1.0,
                        help='Minimum bond length threshold in Angstrom (default: 1.0)')
    parser.add_argument('--base-dir', '-d', default='.',
                        help='Base directory to search for structures.db (default: current directory)')
    args = parser.parse_args()

    formula_dirs = find_formula_dirs_with_db(args.base_dir)
    print(f"Found {len(formula_dirs)} structures.db file(s) under {args.base_dir}")
    print(f"Bond length threshold: {args.threshold} Angstrom")

    if not formula_dirs:
        return

    if 'SLURM_CPUS_PER_TASK' in os.environ:
        num_processes = int(os.environ['SLURM_CPUS_PER_TASK'])
    elif 'PBS_NUM_PPN' in os.environ:
        num_processes = int(os.environ['PBS_NUM_PPN'])
    else:
        num_processes = multiprocessing.cpu_count()
    num_processes = max(1, min(num_processes, len(formula_dirs)))

    print(f"Using {num_processes} processes for parallel processing")

    total_checked = 0
    total_marked = 0
    if num_processes > 1 and len(formula_dirs) > 1:
        with ProcessPoolExecutor(max_workers=num_processes) as ex:
            futures = {
                ex.submit(process_formula_db, fd, args.threshold): fd
                for fd in formula_dirs
            }
            for fut in as_completed(futures):
                fd = futures[fut]
                try:
                    nc, nm = fut.result()
                    total_checked += nc
                    total_marked += nm
                except Exception as exc:  # noqa: BLE001
                    print(f"Error processing {fd}: {exc}")
    else:
        for fd in formula_dirs:
            nc, nm = process_formula_db(fd, args.threshold)
            total_checked += nc
            total_marked += nm

    print(f"\n[bond_check] checked={total_checked}  marked_bond_mis={total_marked}")


if __name__ == '__main__':
    main()
