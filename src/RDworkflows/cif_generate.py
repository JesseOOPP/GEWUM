import logging
from pyxtal import pyxtal
import numpy as np
from pyxtal.msg import Comp_CompatibilityError
import os
import tempfile
from multiprocessing import Pool
import argparse
from ase.io import read, write
from gewum.src.common.cif_db import CifDatabase

DIM_CONFIG = {
    0: {
        'name': 'point group',
        'groups': list(range(2, 57)),
        'default_max_atoms': 36,
        'from_random_kwargs': {'factor': 1},
        'add_vacuum': True,
    },
    1: {
        'name': 'rod group',
        'groups': list(range(2, 76)),
        'default_max_atoms': 36,
        'from_random_kwargs': {'factor': 1},
        'add_vacuum': False,
    },
    2: {
        'name': 'layer group',
        'groups': list(range(2, 81)),
        'default_max_atoms': 36,
        'from_random_kwargs': {'thickness': 3, 'factor': 1},
        'add_vacuum': False,
    },
    3: {
        'name': 'space group',
        'groups': list(range(2, 231)),
        'default_max_atoms': 36,
        'from_random_kwargs': {},
        'add_vacuum': False,
    },
}


def read_input_file(file_path, default_max_atoms):
    inputs = []
    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            data = eval(line.strip())
            if len(data) == 3:
                elements, base_ratio, max_xtals = data
                max_atoms = default_max_atoms
            else:
                elements, base_ratio, max_xtals, max_atoms = data
            inputs.append((elements, base_ratio, int(max_xtals), max_atoms))
    return inputs


def create_folder(elements, base_ratio):
    folder_name = ''.join([f"{elem}{count}" for elem, count in zip(elements, base_ratio)])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    return folder_name


def calculate_possible_ratios(base_ratio, max_atoms):
    possible_ratios = []
    for multiplicity in range(1, max_atoms // sum(base_ratio) + 1):
        atom_counts = tuple(x * multiplicity for x in base_ratio)
        if sum(atom_counts) <= max_atoms:
            possible_ratios.append(atom_counts)
    return possible_ratios


def _pyxtal_to_cif_string(xtal, add_vacuum=False, vacuum=10.0):
    """Convert a pyxtal crystal object to CIF string via temp file (auto-cleaned).

    Uses xtal.to_file() (the only reliable pyxtal output API) with a temp file
    that is read back and immediately deleted.  For 0D, adds vacuum via ASE.

    Returns:
        CIF content string, or None on failure.
    """
    tmp_path = None
    try:
        if add_vacuum:
            with tempfile.NamedTemporaryFile(suffix='.xyz', delete=False) as tmp:
                tmp_path = tmp.name
            xtal.to_file(tmp_path)
            atoms = read(tmp_path, format='xyz')
            atoms.center(vacuum=vacuum)
            with tempfile.NamedTemporaryFile(suffix='.cif', delete=False) as tmp2:
                cif_path = tmp2.name
            write(cif_path, atoms, format='cif')
            with open(cif_path, 'r', encoding='utf-8', errors='replace') as f:
                content = f.read()
            os.remove(cif_path)
            return content
        else:
            with tempfile.NamedTemporaryFile(suffix='.cif', delete=False) as tmp:
                tmp_path = tmp.name
            xtal.to_file(tmp_path)
            with open(tmp_path, 'r', encoding='utf-8', errors='replace') as f:
                return f.read()
    except Exception as e:
        logging.error(f"Failed to convert crystal to CIF string: {e}")
        return None
    finally:
        if tmp_path is not None:
            try:
                os.remove(tmp_path)
            except OSError:
                pass


def generate_crystals_for_group(args):
    elements, base_ratio, max_xtals, max_atoms, group, folder_name, max_attempts, dim, from_random_kwargs, group_name, add_vacuum = args

    possible_ratios = calculate_possible_ratios(base_ratio, max_atoms)

    incompatible_ratios = []
    num_xtals = 0
    skip_group = False
    attempts = 0
    results = [] 

    main_rng = np.random.default_rng()

    while num_xtals < max_xtals and not skip_group and attempts < max_attempts:
        seed = main_rng.integers(0, 2**32)
        sub_rng = np.random.default_rng(seed)
        
        available_ratios = [ratio for ratio in possible_ratios if ratio not in incompatible_ratios]
        if not available_ratios:
            logging.warning(f"No more available ratios for {group_name} {group}. Skipping...")
            break
        
        atom_counts = available_ratios[sub_rng.integers(0, len(available_ratios))]
        
        try:
            xtal = pyxtal()
            xtal.from_random(dim, group, elements, atom_counts, random_state=sub_rng, **from_random_kwargs)

            cif_name = f'xtal_{num_xtals + 1}.cif'
            cif_content = _pyxtal_to_cif_string(xtal, add_vacuum=add_vacuum)

            if cif_content is not None and cif_content.strip():
                results.append((group, cif_name, cif_content))
                num_xtals += 1
                logging.info(f"Generated crystal {num_xtals}/{max_xtals} for {elements} {group_name} {group}")
            else:
                logging.warning(f"Failed to get valid CIF string for crystal in {group_name} {group}")

        except (Comp_CompatibilityError, RuntimeError) as e:
            if "long time to generate structure" in str(e):
                logging.warning(f"{group_name.capitalize()} {group} skipped due to timeout: {str(e)}")
                skip_group = True
                break
            else:
                if atom_counts not in incompatible_ratios:
                    incompatible_ratios.append(atom_counts)
                logging.warning(f"Composition {atom_counts} incompatible for {elements} {group_name} {group}: {str(e)}")
        
        attempts += 1

    if attempts >= max_attempts:
        logging.warning(f"Reached maximum attempts ({max_attempts}) for {elements} {group_name} {group}")
    
    return (elements, group_name, group, num_xtals, max_xtals, results)


def main():
    """Main entry point for crystal generation"""
    parser = argparse.ArgumentParser(description='Generate crystal structures for 0D/1D/2D/3D')
    parser.add_argument('--dim', type=int, required=True, choices=[0, 1, 2, 3],
                       help='Dimension: 0 (point groups), 1 (rod groups), 2 (layer groups), 3 (space groups)')
    parser.add_argument('--max-atoms', type=int, default=None,
                       help='Override maximum atoms for all compositions')
    parser.add_argument('--input-file', type=str, default='cifgen.inp',
                       help='Input file path (default: cifgen.inp)')
    parser.add_argument('--max-attempts', type=int, default=30,
                       help='Maximum attempts per group (default: 30)')
    parser.add_argument('--groups', type=str, default=None,
                       help='Comma-separated space-group numbers to generate '
                            '(e.g. "225,194,139"); default: all groups for the dimension')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    config = DIM_CONFIG[args.dim]
    group_name = config['name']
    groups = config['groups']
    if args.groups is not None:
        requested = [int(g) for g in args.groups.split(',') if g.strip()]
        valid = set(config['groups'])
        groups = [g for g in requested if g in valid]
        skipped = [g for g in requested if g not in valid]
        if skipped:
            logging.warning(f"Ignoring out-of-range groups: {skipped}")
        if not groups:
            logging.error(f"No valid groups in --groups '{args.groups}'")
            return
    default_max_atoms = config['default_max_atoms']
    from_random_kwargs = config['from_random_kwargs']
    add_vacuum = config['add_vacuum']

    inputs = read_input_file(args.input_file, default_max_atoms)

    if args.max_atoms is not None:
        inputs = [(elements, base_ratio, max_xtals, args.max_atoms) 
                 for elements, base_ratio, max_xtals, _ in inputs]

    if 'SLURM_CPUS_PER_TASK' in os.environ:
        max_workers_t = int(os.environ['SLURM_CPUS_PER_TASK'])
    else:
        max_workers_t = os.cpu_count()

    logging.info(f"Generating {args.dim}D crystals ({group_name}s {groups[0]}-{groups[-1]})")
    logging.info(f"Number of CPU cores to use: {max_workers_t}")

    formula_tasks = {}
    seen_folders = set()
    for elements, base_ratio, max_xtals, max_atoms in inputs:
        folder_name = create_folder(elements, base_ratio)
        if folder_name in seen_folders:
            logging.warning(
                f"Duplicate composition '{folder_name}' detected in input; "
                f"skipping the redundant entry to avoid overwriting structures."
            )
            continue
        seen_folders.add(folder_name)
        formula_tasks[folder_name] = []
        for group in groups:
            formula_tasks[folder_name].append(
                (elements, base_ratio, max_xtals, max_atoms, group, folder_name,
                 args.max_attempts, args.dim, from_random_kwargs, group_name, add_vacuum)
            )

    total_tasks = sum(len(t) for t in formula_tasks.values())
    logging.info(f"Total tasks to process: {total_tasks} across {len(formula_tasks)} formulas")
    logging.info(f"Starting parallel processing with {max_workers_t} workers")

    for folder_name, tasks in formula_tasks.items():
        logging.info(f"Generating structures for {folder_name} ({len(tasks)} groups)...")

        all_entries = []
        with Pool(max_workers_t) as pool:
            for worker_result in pool.imap_unordered(generate_crystals_for_group, tasks):
                elements, grp_name, grp, num_xtals, max_xtals, entries = worker_result
                logging.info(f"Completed {elements} {grp_name} {grp}: {num_xtals}/{max_xtals} crystals")
                all_entries.extend(entries)

        abs_folder = os.path.abspath(folder_name)
        os.makedirs(abs_folder, exist_ok=True)

        seen_keys = set()
        unique_entries = []
        for sg, name, content in all_entries:
            key = (sg, name)
            if key in seen_keys:
                logging.warning(
                    f"Duplicate key sg={sg}, name={name} in {folder_name}; "
                    f"the later record will overwrite the earlier one. Skipping."
                )
                continue
            seen_keys.add(key)
            unique_entries.append((sg, name, content))

        try:
            with CifDatabase(abs_folder) as db:
                db.init_db()
                if unique_entries:
                    records = [(sg, name, content, 'initial', None, None, folder_name)
                               for sg, name, content in unique_entries]
                    db.insert_batch(records)
                    logging.info(
                        f"Inserted {len(unique_entries)} structures into "
                        f"{abs_folder}/structures.db"
                    )
                else:
                    logging.warning(
                        f"No structures generated for {folder_name}; "
                        f"created empty structures.db"
                    )
        except Exception:
            logging.exception(
                f"Failed to write SQLite DB for {folder_name}; "
                f"continuing with remaining formulas."
            )

    logging.info(f"All {args.dim}D crystals generated and stored in SQLite databases.")


if __name__ == "__main__":
    main()
