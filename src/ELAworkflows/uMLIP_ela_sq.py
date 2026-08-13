import os
import glob
import argparse
import numpy as np
from ase.io import read, write
from ase.optimize import BFGS
from ase.filters import UnitCellFilter
from ase.units import GPa
from mattersim.forcefield import MatterSimCalculator

MODEL_PATH = "MatterSim-v1.0.0-1M.pth"


def relax(atoms, calc, pressure=0.0, maxstep=0.1, eps=None, max_step=None):
    atoms.calc = calc
    mask = [False, False, False, False, False, False]
    ucf = UnitCellFilter(atoms, scalar_pressure=pressure * GPa, mask=mask, constant_volume=False)
    dyn = BFGS(ucf, maxstep=maxstep)
    dyn.run(fmax=eps, steps=max_step)
    return atoms


def _write_outcar(outcar_path, device, energy, max_force, fmax):
    """Write an OUTCAR that mimics the VASP fields vaspkit -task 201 parses."""
    with open(outcar_path, 'w') as f:
        f.write(f"Using MatterSim on {device} device\n")
        if max_force <= fmax:
            f.write(f"  free  energy   TOTEN  = {energy:.6f} eV\n")
            f.write("                 Voluntary context switches:\n")
            f.write("Optimized structure saved to CONTCAR\n")
        else:
            f.write(f"Warning: Optimization did not fully converge! Max force: {max_force:.6f} eV/A\n")


def process_one(poscar_path, contcar_path, outcar_path, calc, device, fmax, max_steps):
    """Relax a single strained structure and emit OUTCAR + CONTCAR in place."""
    atoms = read(poscar_path)
    relaxed_atoms = relax(
        atoms=atoms,
        calc=calc,
        pressure=0.0,
        maxstep=0.1,
        eps=fmax,
        max_step=max_steps,
    )

    atom_force = relaxed_atoms.get_forces()
    max_force = float(np.max(np.abs(atom_force)))
    energy = float(relaxed_atoms.get_potential_energy())

    _write_outcar(outcar_path, device, energy, max_force, fmax)
    write(contcar_path, relaxed_atoms, vasp5=True, direct=True)
    return max_force <= fmax


def run_batch(base_dir, calc, device, fmax, max_steps):
    """Load the model once (by the caller) and process every strain_* directory."""
    strain_dirs = sorted(
        d for d in glob.glob(os.path.join(base_dir, '**', 'strain_*'), recursive=True)
        if os.path.isdir(d)
    )

    if not strain_dirs:
        print(f"No strain_* directories found under {base_dir}")
        return

    print(f"Found {len(strain_dirs)} strain directories under {base_dir}")
    for strain_dir in strain_dirs:
        poscar_path = os.path.join(strain_dir, 'POSCAR')
        if not os.path.exists(poscar_path):
            print(f"Skip {strain_dir}: no POSCAR")
            continue
        try:
            converged = process_one(
                poscar_path,
                os.path.join(strain_dir, 'CONTCAR'),
                os.path.join(strain_dir, 'OUTCAR'),
                calc,
                device,
                fmax,
                max_steps,
            )
            status = "converged" if converged else "NOT converged"
            print(f"Done {strain_dir}: {status}")
        except Exception as exc:
            print(f"Error {strain_dir}: {exc}")


def main():
    parser = argparse.ArgumentParser(description='GEWUM ELA Single-point Calculation')
    parser.add_argument('--input', '-i', default='POSCAR', help='Input POSCAR file (single mode)')
    parser.add_argument('--output', '-o', default='CONTCAR', help='Output CONTCAR file (single mode)')
    parser.add_argument('--batch', default=None,
                        help='Base directory: relax every strain_* subdirectory with a single model load')
    parser.add_argument('--fmax', type=float, default=0.008, help='Force convergence threshold')
    parser.add_argument('--max-steps', type=int, default=300, help='Maximum optimization steps')
    parser.add_argument('--device', default='cpu', help='Device for calculation')
    args = parser.parse_args()

    print(f"Using MatterSim on {args.device} device")
    calc = MatterSimCalculator(load_path=MODEL_PATH, device=args.device)

    if args.batch:
        run_batch(args.batch, calc, args.device, args.fmax, args.max_steps)
    else:
        converged = process_one(
            args.input, args.output, 'OUTCAR',
            calc, args.device, args.fmax, args.max_steps,
        )
        if converged:
            print(f"Optimized structure saved to {args.output}")
        else:
            print("Warning: Optimization did not fully converge!")


if __name__ == "__main__":
    main()
