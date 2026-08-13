"""
GEWUM SQ AM Module
Stochastic Quenching for MD-free amorphous structure generation.

Generates amorphous configurations by placing atoms randomly at the target
density (subject to a minimum pairwise-distance constraint) and directly
relaxing them to the nearest local minimum with MatterSim, without any
melt-quench molecular dynamics. Repeating over independent random seeds
yields an ensemble of glassy structures.

The output artifacts are a directory of relaxed CIFs plus am_ml_energy.csv
(one row per configuration: frame name, potential energy and per-atom energy).

Reference: Holmstrom et al., Phys. Rev. B 79, 144201 (2009).
"""
import os
import csv
import sys
import logging
import argparse
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.build import make_supercell
from ase.data import covalent_radii
from ase.optimize import FIRE, LBFGS
from ase.filters import UnitCellFilter
from mattersim.forcefield import MatterSimCalculator

# Reuse the PT-module perturbation operators (read-only cross-module import).
# These provide rattle + optional lattice strain + optional atom rotation plus
# min-distance/uniqueness checks, identical to the PT workflow's perturbation.
from gewum.src.cifgen_input.mutate import (
    generate_modified_structure,
    calculate_structure_fingerprint,
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Default MatterSim model (matches the AM featurization pipeline)
DEFAULT_MODEL = "mattersim-v1.0.0-5M.pth"

# CSV schema for the per-frame energy table
CSV_FIELDS = ['frame', 'potential_energy_eV', 'energy_per_atom_eV',
              'temperature_K', 'stage', 'step', 'converged', 'max_force']


def create_supercell(atoms, supercell_size):
    """Expand a unit cell to a supercell."""
    P = np.diag(supercell_size)
    supercell = make_supercell(atoms, P)
    logging.info(f"Supercell: {supercell_size[0]}x{supercell_size[1]}x{supercell_size[2]}")
    logging.info(f"Atoms: {len(atoms)} -> {len(supercell)}")
    return supercell


def _min_image_distance(cart_delta, cell, inv_cell):
    """Minimum-image distance for a cartesian displacement in a periodic cell."""
    frac = cart_delta @ inv_cell
    frac -= np.round(frac)
    dmic = frac @ cell
    return float(np.linalg.norm(dmic))


def generate_random_config(template, d_min_scale=0.7, max_tries=20000, rng=None):
    """
    Build a random atomic configuration with the same composition and cell
    as `template`, placing atoms one at a time and rejecting any placement
    that violates the minimum pairwise-distance constraint.

    Args:
        template: ASE Atoms providing composition, cell and density
        d_min_scale: minimum pairwise distance = d_min_scale*(r_cov_i + r_cov_j)
        max_tries: max rejection-sampling attempts per atom before giving up
        rng: numpy random Generator (for reproducibility)

    Returns:
        ASE Atoms with randomized fractional positions, or None if the box is
        too dense to satisfy the distance constraint.
    """
    if rng is None:
        rng = np.random.default_rng()

    numbers = template.get_atomic_numbers()
    cell = np.array(template.get_cell())
    inv_cell = np.linalg.inv(cell)
    radii = np.array([covalent_radii[z] for z in numbers])

    placed_cart = []
    placed_idx = []

    # Place atoms in descending-radius order so large atoms find room first
    order = np.argsort(-radii)

    for i in order:
        r_i = radii[i]
        success = False
        for _ in range(max_tries):
            frac = rng.random(3)
            cart = frac @ cell
            ok = True
            for j, cart_j in zip(placed_idx, placed_cart):
                d_min = d_min_scale * (r_i + radii[j])
                if _min_image_distance(cart - cart_j, cell, inv_cell) < d_min:
                    ok = False
                    break
            if ok:
                placed_cart.append(cart)
                placed_idx.append(i)
                success = True
                break
        if not success:
            logging.error(
                f"Failed to place atom {i} (Z={numbers[i]}) after {max_tries} tries; "
                f"box likely too dense. Try a larger supercell or smaller --d-min-scale.")
            return None

    # Reassemble positions back into original atom order
    positions = np.zeros((len(numbers), 3))
    for idx, cart in zip(placed_idx, placed_cart):
        positions[idx] = cart

    atoms = Atoms(numbers=numbers, positions=positions, cell=cell, pbc=True)
    return atoms


def relax_config(atoms, calc, variable_cell=True, fmax=0.05,
                 steps_fire=300, steps_final=200):
    """
    Relax a random configuration to its nearest local minimum.

    A robust two-stage scheme is used: FIRE first (tolerant of the large,
    out-of-distribution forces from random starts), then LBFGS refinement.
    When variable_cell is True the cell is relaxed at zero pressure so the
    density self-adjusts.

    Returns:
        (relaxed_atoms, energy, converged, max_force) or None on failure.
    """
    atoms.calc = calc
    try:
        target = UnitCellFilter(atoms) if variable_cell else atoms

        dyn = FIRE(target, logfile=None)
        dyn.run(fmax=fmax, steps=steps_fire)

        if steps_final > 0:
            dyn2 = LBFGS(target, logfile=None)
            dyn2.run(fmax=fmax, steps=steps_final)

        energy = float(atoms.get_potential_energy())
        if not np.isfinite(energy):
            logging.error("Relaxation produced a non-finite energy; skipping.")
            return None

        forces = atoms.get_forces()
        max_force = float(np.max(np.linalg.norm(forces, axis=1)))
        converged = max_force <= fmax * 1.5
        return atoms, energy, converged, max_force
    except Exception as e:
        logging.error(f"Relaxation failed: {e}")
        return None
    finally:
        atoms.calc = None


def run_sq(cif_file, output_dir="sq_output", n_configs=20,
           supercell_size=(2, 2, 2), variable_cell=True, fmax=0.05,
           d_min_scale=0.7, steps_fire=300, steps_final=200,
           model_path=DEFAULT_MODEL, device="cpu", seed_base=0,
           init_mode="random", perturb_params=None, save_initial=False):
    """
    Run stochastic quenching to generate an amorphous ensemble.

    The MatterSim model is loaded once and reused across all configurations.

    init_mode selects how each initial configuration is built:
        "random"  - fully random packing (unbiased, no crystal memory)
        "perturb" - perturb the crystal template via the PT-module operators
                    (rattle + optional lattice strain + optional rotation),
                    inheriting the template's structural motif
    perturb_params is the dict consumed by generate_modified_structure and is
    only required when init_mode == "perturb".

    Returns:
        dict summary, or None on fatal failure.
    """
    try:
        template = read(cif_file)
        base_name = os.path.basename(cif_file).replace('.cif', '')
        logging.info(f"Starting stochastic quenching for: {base_name}")
        logging.info(f"Original structure: {len(template)} atoms")

        if any(s > 1 for s in supercell_size):
            template = create_supercell(template, supercell_size)

        # Load the MatterSim model once and reuse across all configurations
        if model_path:
            calc = MatterSimCalculator(load_path=model_path, device=device)
        else:
            calc = MatterSimCalculator(device=device)

        os.makedirs(output_dir, exist_ok=True)
        frame_dir = os.path.join(output_dir, "am_frames")
        os.makedirs(frame_dir, exist_ok=True)
        # Optional: keep the pre-relaxation initial configs in init_frames/.
        # An initial is named <base>_sqNNNN_init.cif and pairs with the relaxed
        # frame <base>_sqNNNN.cif in am_frames/ (same NNNN index). It is written
        # only after a successful relaxation, so the two dirs stay one-to-one.
        init_dir = os.path.join(output_dir, "init_frames")
        if save_initial:
            os.makedirs(init_dir, exist_ok=True)
        energy_csv = os.path.join(output_dir, "am_ml_energy.csv")

        n_atoms = len(template)
        energy_records = []
        n_success = 0
        n_fail = 0

        # Origin tag written to the CSV 'stage' column for traceability.
        stage_tag = 'sq_pert' if init_mode == 'perturb' else 'sq'

        # Perturbation route reuses the PT operators; keep a running list of
        # fingerprints (seeded with the pristine template) so relaxed inputs
        # are not near-duplicates of the crystal or of each other.
        tol = 0.1
        existing_fingerprints = []
        if init_mode == 'perturb':
            if perturb_params is None:
                logging.error("init_mode='perturb' requires perturb_params.")
                return None
            tol = perturb_params.get('similarity_tolerance', 0.1)
            existing_fingerprints.append(
                calculate_structure_fingerprint(template, tol))

        for idx in range(n_configs):
            seed = seed_base + idx

            if init_mode == 'perturb':
                # mutate operators use the legacy global RNG; seed per-config
                np.random.seed(seed % (2**32 - 1))
                config = generate_modified_structure(
                    template, perturb_params, existing_fingerprints, tol)
                if config is not None:
                    existing_fingerprints.append(
                        calculate_structure_fingerprint(config, tol))
            else:
                rng = np.random.default_rng(seed)
                config = generate_random_config(
                    template, d_min_scale=d_min_scale, rng=rng)

            if config is None:
                n_fail += 1
                continue

            frame_name = f"{base_name}_sq{idx:04d}.cif"
            # Snapshot the TRUE initial before relaxation (relax_config mutates
            # the atoms in place). It is written only after a successful
            # relaxation so init_frames and am_frames stay strictly one-to-one.
            init_config = config.copy() if save_initial else None

            result = relax_config(
                config, calc, variable_cell=variable_cell, fmax=fmax,
                steps_fire=steps_fire, steps_final=steps_final)
            if result is None:
                n_fail += 1
                continue

            relaxed, energy, converged, max_force = result
            write(os.path.join(frame_dir, frame_name), relaxed)
            if save_initial:
                init_name = f"{base_name}_sq{idx:04d}_init.cif"
                write(os.path.join(init_dir, init_name), init_config)

            energy_records.append({
                'frame': frame_name,
                'potential_energy_eV': f"{energy:.6f}",
                'energy_per_atom_eV': f"{energy / len(relaxed):.6f}",
                'temperature_K': 0,
                'stage': stage_tag,
                'step': idx,
                'converged': converged,
                'max_force': f"{max_force:.6f}",
            })
            n_success += 1
            logging.info(
                f"  config {idx:04d}: E={energy:.4f} eV "
                f"({energy / len(relaxed):.4f} eV/atom), "
                f"converged={converged}")

        if n_success == 0:
            logging.error("All stochastic-quenching configurations failed.")
            return None

        with open(energy_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
            writer.writeheader()
            writer.writerows(energy_records)

        logging.info(f"Stochastic quenching completed: {base_name}")
        logging.info(f"  Success: {n_success}/{n_configs} (failed: {n_fail})")
        logging.info(f"  Frames dir: {frame_dir}")
        logging.info(f"  Energy CSV: {energy_csv}")

        return {
            'frame_count': n_success,
            'fail_count': n_fail,
            'output_dir': output_dir,
            'frame_dir': frame_dir,
            'init_dir': init_dir if save_initial else None,
            'energy_csv': energy_csv,
            'n_atoms': n_atoms,
        }

    except Exception as e:
        logging.error(f"Stochastic quenching failed for {cif_file}: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return None


def main():
    parser = argparse.ArgumentParser(
        description='GEWUM SQ - Stochastic Quenching for MD-free amorphous generation',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('cif_file', help='Input CIF file path (crystal template)')
    parser.add_argument('--output-dir', '-o', default='sq_output',
                        help='Output directory (default: sq_output)')
    parser.add_argument('--n-configs', '-n', type=int, default=20,
                        help='Number of random configurations (default: 20)')
    parser.add_argument('--supercell', type=int, nargs=3, default=[2, 2, 2],
                        metavar=('NX', 'NY', 'NZ'),
                        help='Supercell expansion for target atom count/density (default: 2 2 2)')
    parser.add_argument('--fixed-cell', action='store_true',
                        help='Relax ionic positions only at fixed cell '
                             '(default: variable cell at zero pressure)')
    parser.add_argument('--fmax', type=float, default=0.05,
                        help='Force convergence threshold in eV/Angstrom (default: 0.05)')
    parser.add_argument('--d-min-scale', type=float, default=0.7,
                        help='Min pairwise distance = scale*(r_cov_i+r_cov_j) (default: 0.7)')
    parser.add_argument('--steps-fire', type=int, default=300,
                        help='FIRE coarse relaxation steps (default: 300)')
    parser.add_argument('--steps-final', type=int, default=200,
                        help='LBFGS refinement steps (default: 200)')
    parser.add_argument('--model', type=str, default=DEFAULT_MODEL,
                        help=f'Path to MatterSim model file (default: {DEFAULT_MODEL})')
    parser.add_argument('--device', type=str, default='cpu',
                        choices=['cpu', 'cuda'],
                        help='Computation device (default: cpu)')
    parser.add_argument('--seed', type=int, default=0,
                        help='Base random seed; config i uses seed+i (default: 0)')
    parser.add_argument('--init-mode', type=str, default='random',
                        choices=['random', 'perturb'],
                        help='Initial-configuration source: fully random packing '
                             'or perturbation of the crystal template (default: random)')
    parser.add_argument('--rattle-stdev-min', type=float, default=0.4,
                        help='Perturb route: min rattle stdev in Angstrom (default: 0.4)')
    parser.add_argument('--rattle-stdev-max', type=float, default=1.0,
                        help='Perturb route: max rattle stdev in Angstrom (default: 1.0)')
    parser.add_argument('--strain-prob', type=float, default=0.5,
                        help='Perturb route: probability of applying lattice strain (default: 0.5)')
    parser.add_argument('--max-strain', type=float, default=0.05,
                        help='Perturb route: max lattice strain magnitude (default: 0.05)')
    parser.add_argument('--rotation-prob', type=float, default=0.3,
                        help='Perturb route: probability of applying atom rotation (default: 0.3)')
    parser.add_argument('--rotation-per-atom', type=float, default=0.15,
                        help='Perturb route: per-atom rotation probability (default: 0.15)')
    parser.add_argument('--max-rotation-angle', type=float, default=60.0,
                        help='Perturb route: max rotation angle in degrees (default: 60.0)')
    parser.add_argument('--perturb-min-dist', type=float, default=1.0,
                        help='Perturb route: min pairwise distance after perturbation in Angstrom (default: 1.0)')
    parser.add_argument('--save-initial', action='store_true',
                        help='Also save the pre-relaxation initial configs to '
                             'init_frames/ (same filename as the relaxed frame)')
    parser.add_argument('--from-config', type=str, default=None,
                        help='Load SQ parameters from am_config.yaml (sq: section)')

    args = parser.parse_args()

    if args.from_config:
        try:
            import yaml
            with open(args.from_config, 'r') as f:
                cfg = yaml.safe_load(f) or {}
            sq = cfg.get('sq', {}) or {}
            args.n_configs = sq.get('n_configs', args.n_configs)
            args.supercell = sq.get('supercell', args.supercell)
            args.fmax = sq.get('fmax', args.fmax)
            args.d_min_scale = sq.get('d_min_scale', args.d_min_scale)
            args.steps_fire = sq.get('steps_fire', args.steps_fire)
            args.steps_final = sq.get('steps_final', args.steps_final)
            args.model = sq.get('model_path', cfg.get('model_path', args.model))
            if 'variable_cell' in sq:
                args.fixed_cell = not bool(sq['variable_cell'])
            logging.info(f"Loaded SQ parameters from {args.from_config}")
        except Exception as e:
            logging.warning(f"Failed to load config: {e}, using CLI defaults")

    perturb_params = {
        'stdev_min': args.rattle_stdev_min,
        'stdev_max': args.rattle_stdev_max,
        'lattice_strain_probability': args.strain_prob,
        'max_strain': args.max_strain,
        'atom_rotation_probability': args.rotation_prob,
        'rotation_probability_per_atom': args.rotation_per_atom,
        'max_rotation_angle': args.max_rotation_angle,
        'min_distance': args.perturb_min_dist,
        'similarity_tolerance': 0.1,
    }

    result = run_sq(
        cif_file=args.cif_file,
        output_dir=args.output_dir,
        n_configs=args.n_configs,
        supercell_size=tuple(args.supercell),
        variable_cell=not args.fixed_cell,
        fmax=args.fmax,
        d_min_scale=args.d_min_scale,
        steps_fire=args.steps_fire,
        steps_final=args.steps_final,
        model_path=args.model,
        device=args.device,
        seed_base=args.seed,
        init_mode=args.init_mode,
        perturb_params=perturb_params,
        save_initial=args.save_initial,
    )
    sys.exit(0 if result else 1)


if __name__ == "__main__":
    main()
