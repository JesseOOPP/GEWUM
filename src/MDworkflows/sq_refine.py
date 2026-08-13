"""
GEWUM SQ AM Refinement

Final-stage precise relaxation for the maxmin-selected amorphous subset.

The advanced SQ pipeline coarse-relaxes a large pool cheaply, screens out
recrystallized frames, then reduces the clean amorphous pool to a diverse
subset. This module spends the expensive LBFGS refinement only on that subset,
reusing the exact relaxation physics of sq_am (variable-cell zero-pressure
FIRE + LBFGS) so there is a single source of truth for the relaxation.

Input is a directory of coarse-relaxed CIF frames (top level only; the
crystalline/ and remove/ sub-directories produced by the earlier stages are
ignored). Output mirrors sq_am: <output_dir>/am_frames/<name>.cif plus
<output_dir>/am_ml_energy.csv with stage='sq_final'.
"""
import os
import csv
import sys
import glob
import logging
import argparse
from ase.io import read, write
from mattersim.forcefield import MatterSimCalculator

# Reuse the SQ relaxation physics and CSV schema (single source of truth).
from gewum.src.MDworkflows.sq_am import relax_config, CSV_FIELDS, DEFAULT_MODEL

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

STAGE_TAG = 'sq_final'


def refine_directory(input_dir, output_dir="am_final", variable_cell=True,
                     fmax=0.05, steps_fire=50, steps_final=200,
                     model_path=DEFAULT_MODEL, device="cpu",
                     shard_index=0, shard_count=1,
                     energy_csv_name="am_ml_energy.csv"):
    """Refine every top-level CIF in input_dir with a final LBFGS relaxation.

    The MatterSim model is loaded once and reused across all frames. Frames
    whose output already exists are skipped (resume-safe).

    For parallel refinement the selected pool is sharded across processes:
    with shard_count > 1 this process only handles files whose sorted index
    satisfies (index %% shard_count == shard_index). All shards write frames
    into the SAME output_dir/am_frames (distinct file names, safe), but each
    writes its own energy_csv_name so the caller can merge them afterwards.

    Returns:
        dict summary, or None on fatal failure.
    """
    if not os.path.isdir(input_dir):
        logging.error(f"Input directory not found: {input_dir}")
        return None

    all_cifs = sorted(glob.glob(os.path.join(input_dir, '*.cif')))
    if not all_cifs:
        logging.error(f"No CIF frames to refine in {input_dir}")
        return None

    if shard_count > 1:
        cif_files = [f for i, f in enumerate(all_cifs)
                     if i % shard_count == shard_index]
        logging.info(
            f"Shard {shard_index}/{shard_count}: "
            f"{len(cif_files)}/{len(all_cifs)} frames")
    else:
        cif_files = all_cifs

    os.makedirs(output_dir, exist_ok=True)
    frame_dir = os.path.join(output_dir, "am_frames")
    os.makedirs(frame_dir, exist_ok=True)
    energy_csv = os.path.join(output_dir, energy_csv_name)

    if not cif_files:
        # Empty shard (more workers than frames) is not an error.
        logging.info("This shard has no frames to refine; nothing to do.")
        return {
            'frame_count': 0, 'skip_count': 0, 'fail_count': 0,
            'output_dir': output_dir, 'frame_dir': frame_dir,
            'energy_csv': energy_csv,
        }

    try:
        if model_path:
            calc = MatterSimCalculator(load_path=model_path, device=device)
        else:
            calc = MatterSimCalculator(device=device)
    except Exception as e:
        logging.error(f"Failed to load MatterSim model: {e}")
        return None

    energy_records = []
    n_success = 0
    n_fail = 0
    n_skip = 0

    for step, path in enumerate(cif_files):
        frame_name = os.path.basename(path)
        out_path = os.path.join(frame_dir, frame_name)
        if os.path.exists(out_path):
            n_skip += 1
            continue

        try:
            atoms = read(path)
        except Exception as e:
            logging.warning(f"Skipping unreadable frame {frame_name}: {e}")
            n_fail += 1
            continue

        result = relax_config(
            atoms, calc, variable_cell=variable_cell, fmax=fmax,
            steps_fire=steps_fire, steps_final=steps_final)
        if result is None:
            n_fail += 1
            continue

        relaxed, energy, converged, max_force = result
        write(out_path, relaxed)
        energy_records.append({
            'frame': frame_name,
            'potential_energy_eV': f"{energy:.6f}",
            'energy_per_atom_eV': f"{energy / len(relaxed):.6f}",
            'temperature_K': 0,
            'stage': STAGE_TAG,
            'step': step,
            'converged': converged,
            'max_force': f"{max_force:.6f}",
        })
        n_success += 1
        logging.info(
            f"  refine {frame_name}: E={energy:.4f} eV "
            f"({energy / len(relaxed):.4f} eV/atom), converged={converged}")

    if n_success == 0 and n_skip == 0:
        logging.error("All refinement relaxations failed.")
        return None

    # Merge with any previous records (resume) so the CSV stays complete.
    if os.path.exists(energy_csv):
        try:
            with open(energy_csv, 'r', newline='') as f:
                prev = list(csv.DictReader(f))
            done = {r['frame'] for r in energy_records}
            energy_records = [r for r in prev if r['frame'] not in done] + energy_records
        except Exception as e:
            logging.warning(f"Could not merge existing CSV: {e}")

    with open(energy_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(energy_records)

    logging.info("SQ refinement completed.")
    logging.info(f"  Success: {n_success}, skipped(resume): {n_skip}, failed: {n_fail}")
    logging.info(f"  Frames dir: {frame_dir}")
    logging.info(f"  Energy CSV: {energy_csv}")

    return {
        'frame_count': n_success,
        'skip_count': n_skip,
        'fail_count': n_fail,
        'output_dir': output_dir,
        'frame_dir': frame_dir,
        'energy_csv': energy_csv,
    }


def main():
    parser = argparse.ArgumentParser(
        description='GEWUM SQ - final LBFGS refinement of the selected subset',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('input_dir',
                        help='Directory of coarse-relaxed CIF frames (top level)')
    parser.add_argument('--output-dir', '-o', default='am_final',
                        help='Output directory (default: am_final)')
    parser.add_argument('--fixed-cell', action='store_true',
                        help='Relax ionic positions only at fixed cell '
                             '(default: variable cell at zero pressure)')
    parser.add_argument('--fmax', type=float, default=0.05,
                        help='Force convergence threshold in eV/Angstrom (default: 0.05)')
    parser.add_argument('--steps-fire', type=int, default=50,
                        help='FIRE steps before LBFGS (default: 50)')
    parser.add_argument('--steps-final', type=int, default=200,
                        help='LBFGS refinement steps (default: 200)')
    parser.add_argument('--model', type=str, default=DEFAULT_MODEL,
                        help=f'Path to MatterSim model file (default: {DEFAULT_MODEL})')
    parser.add_argument('--device', type=str, default='cpu',
                        choices=['cpu', 'cuda'],
                        help='Computation device (default: cpu)')
    parser.add_argument('--shard-index', type=int, default=0,
                        help='Index of this shard (0-based) for parallel refine')
    parser.add_argument('--shard-count', type=int, default=1,
                        help='Total number of parallel shards (default: 1)')
    parser.add_argument('--energy-csv-name', type=str, default='am_ml_energy.csv',
                        help='Energy CSV file name inside output-dir '
                             '(per-shard name for parallel refine)')

    args = parser.parse_args()

    result = refine_directory(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        variable_cell=not args.fixed_cell,
        fmax=args.fmax,
        steps_fire=args.steps_fire,
        steps_final=args.steps_final,
        model_path=args.model,
        device=args.device,
        shard_index=args.shard_index,
        shard_count=args.shard_count,
        energy_csv_name=args.energy_csv_name,
    )
    sys.exit(0 if result else 1)


if __name__ == "__main__":
    main()
