"""
GEWUM Unified Structure Relaxation Module
Supports RD (Random Design), PT (Perturbation), and HP (High Pressure) workflows
"""
import os
import logging
import csv
import gc
import argparse
import tempfile
from ase.io import read, write
from ase.optimize import BFGS, LBFGS
from ase.filters import UnitCellFilter
from mattersim.forcefield import MatterSimCalculator
#from upet.calculator import UPETCalculator
#from deepmd.calculator import DP as DPCalculator

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

GPA_TO_EV_PER_ANG3 = 1 / 160.21766208  # 1 GPa = 1/160.2177 eV/AA3

CSV_HEADER_DEFAULT = ['CIF_Base_Name', 'Total_Energy_eV',
                      'Energy_per_Atom_eV', 'Relaxed_CIF_Path']
CSV_HEADER_HP = ['CIF_Base_Name', 'Total_Energy_eV', 'Energy_per_Atom_eV',
                 'Final_Pressure_GPa', 'Enthalpy_per_Atom_eV',
                 'Corrected_Enthalpy_per_Atom_eV', 'Relaxed_CIF_Path']


def _run_optimization(atoms, mode, fmax, max_steps, pressure):
    """Core in-memory optimization. Returns (atoms, energy, energy_per_atom, hp_data).

    Caller owns the lifetime of `atoms` and its calculator.
    """
    hp_data = None
    dyn = None
    ucf = None
    try:
        if mode == 1:
            dyn = BFGS(atoms)
            dyn.run(fmax=fmax, steps=max_steps)
        elif mode == 2:
            mask = [True, True, True, False, False, False]
            ucf = UnitCellFilter(atoms, mask=mask)
            dyn = LBFGS(ucf)
            dyn.run(fmax=fmax, steps=max_steps)
        elif mode == 3:
            target_pressure = pressure * GPA_TO_EV_PER_ANG3
            ucf = UnitCellFilter(atoms, scalar_pressure=target_pressure)
            dyn = LBFGS(ucf)
            dyn.run(fmax=fmax, steps=max_steps)
            atoms = ucf.atoms
        else:
            raise ValueError(f"Unknown mode: {mode}. Use 1 (atoms only), 2 (atoms + cell), or 3 (high pressure).")

        optimized_energy = float(atoms.get_potential_energy())
        optimized_energy_per_atom = float(optimized_energy / len(atoms))

        if mode == 3 and pressure != 0.0:
            stress_voigt = atoms.calc.get_stress()
            stress_voigt_GPa = stress_voigt * 160.21766208
            hydrostatic_pressure = float(-(stress_voigt_GPa[0] + stress_voigt_GPa[1] + stress_voigt_GPa[2]) / 3)

            volume = float(atoms.get_volume())
            PV_term = hydrostatic_pressure * volume / 160.21766208
            enthalpy = optimized_energy + PV_term
            enthalpy_per_atom = float(enthalpy / len(atoms))

            dE_dPV = float((pressure - hydrostatic_pressure) * volume / 160.21766208 / len(atoms))
            corrected_enthalpy_per_atom = float(enthalpy_per_atom + dE_dPV)

            hp_data = {
                'final_pressure': hydrostatic_pressure,
                'enthalpy_per_atom': enthalpy_per_atom,
                'corrected_enthalpy_per_atom': corrected_enthalpy_per_atom,
                'dE_dPV': dE_dPV,
            }
            logging.info(f"Target pressure: {pressure} GPa, Final pressure: {hydrostatic_pressure:.2f} GPa")
            logging.info(f"Corrected enthalpy: {corrected_enthalpy_per_atom:.6f} eV/atom")

        return atoms, optimized_energy, optimized_energy_per_atom, hp_data
    finally:
        if dyn is not None:
            del dyn
        if ucf is not None:
            del ucf


def _atoms_to_cif_text(atoms):
    """Serialize ASE Atoms to a CIF text string via a temp file (cross-platform safe)."""
    fd, tmp_path = tempfile.mkstemp(suffix='.cif')
    os.close(fd)
    try:
        write(tmp_path, atoms)
        with open(tmp_path, 'r', encoding='utf-8') as fh:
            return fh.read()
    finally:
        try:
            os.remove(tmp_path)
        except OSError:
            pass


def _cif_text_to_atoms(cif_text):
    """Load ASE Atoms from a CIF text string via a temp file."""
    fd, tmp_path = tempfile.mkstemp(suffix='.cif')
    os.close(fd)
    try:
        with open(tmp_path, 'w', encoding='utf-8') as fh:
            fh.write(cif_text)
        return read(tmp_path)
    finally:
        try:
            os.remove(tmp_path)
        except OSError:
            pass


def append_energy_csv(csv_path, mode, base_name, energy, energy_per_atom,
                      hp_data, outfile):
    """Append a single relaxation result row to the per-formula CSV.

    Writes the appropriate header on first append. The caller is responsible
    for serializing concurrent writes (we recommend a single-process commit step).
    """
    file_exists = os.path.isfile(csv_path)
    os.makedirs(os.path.dirname(csv_path) or '.', exist_ok=True)
    with open(csv_path, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        if not file_exists:
            writer.writerow(CSV_HEADER_HP if mode == 3 else CSV_HEADER_DEFAULT)
        if mode == 3 and hp_data:
            writer.writerow([base_name, energy, energy_per_atom,
                             hp_data['final_pressure'], hp_data['enthalpy_per_atom'],
                             hp_data['corrected_enthalpy_per_atom'], outfile])
        else:
            writer.writerow([base_name, energy, energy_per_atom, outfile])


def optimize_from_content(cif_text, base_name, mode=2, fmax=0.05,
                          max_steps=200, pressure=0.0):
    """In-memory relaxation entry point used by the DB-driven worker.

    Args:
        cif_text: Source CIF as a string (typically loaded from structures.db)
        base_name: Logical name for logging (no .cif extension)
        mode/fmax/max_steps/pressure: Same as optimize_structure

    Returns:
        (relaxed_cif_text, energy, energy_per_atom, hp_data)
    """
    device = "cpu"
    if mode == 3:
        logging.info(f"[{base_name}] uMLIP {device} mode={mode} (HP) P={pressure}GPa fmax={fmax}")
    else:
        logging.info(f"[{base_name}] uMLIP {device} mode={mode} fmax={fmax}")

    atoms = None
    calculator = None
    try:
        atoms = _cif_text_to_atoms(cif_text)
        calculator = MatterSimCalculator(device=device)
        atoms.calc = calculator

        atoms, energy, energy_per_atom, hp_data = _run_optimization(
            atoms, mode, fmax, max_steps, pressure
        )

        relaxed_text = _atoms_to_cif_text(atoms)
        return relaxed_text, energy, energy_per_atom, hp_data
    finally:
        if atoms is not None:
            atoms.calc = None
        if calculator is not None:
            del calculator
        gc.collect()


def optimize_structure(cif_file, relaxed_dir, mode=1, fmax=0.05, max_steps=200, pressure=0.0):
    """Legacy disk-based optimization (kept for backward compatibility).

    Reads cif_file from disk, writes {base_name}_relaxed.cif into relaxed_dir,
    and appends to ../energy_results.csv. Returns the CSV path.
    """
    device = "cpu"
    if mode == 3:
        logging.info(f"Running uMLIP on {device}, mode={mode} (HP), pressure={pressure} GPa, fmax={fmax}")
    else:
        logging.info(f"Running uMLIP on {device}, mode={mode}, fmax={fmax}")

    atoms = None
    calculator = None
    try:
        atoms = read(cif_file)
        calculator = MatterSimCalculator(device=device)
        atoms.calc = calculator

        atoms, optimized_energy, optimized_energy_per_atom, hp_data = _run_optimization(
            atoms, mode, fmax, max_steps, pressure
        )

        base_name = os.path.basename(cif_file).replace('.cif', '')
        outfile = os.path.join(relaxed_dir, f"{base_name}_relaxed.cif")
        write(outfile, atoms)

        csv_path = os.path.join(os.path.dirname(relaxed_dir), "energy_results.csv")
        append_energy_csv(csv_path, mode, base_name,
                          optimized_energy, optimized_energy_per_atom,
                          hp_data, outfile)

        logging.info(f"Optimized structure saved to {outfile}")
        logging.info(f"Results appended to CSV: {csv_path}")
        return csv_path
    finally:
        if atoms is not None:
            atoms.calc = None
        if calculator is not None:
            del calculator
        gc.collect()


def main():
    """Command line interface for structure relaxation"""
    parser = argparse.ArgumentParser(
        description='GEWUM Structure Relaxation using MatterSim',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('cif_file', help='Path to input CIF file')
    parser.add_argument('relaxed_dir', help='Directory to save relaxed structure')
    parser.add_argument('--mode', type=int, default=1, choices=[1, 2, 3],
                        help='Optimization mode:\n'
                             '  1 = atoms only (BFGS, for PT workflow)\n'
                             '  2 = atoms + cell (LBFGS+UnitCellFilter, for RD workflow)\n'
                             '  3 = high pressure (LBFGS+UnitCellFilter, for HP workflow)')
    parser.add_argument('--fmax', type=float, default=0.05,
                        help='Force convergence threshold in eV/Angstrom (default: 0.05)')
    parser.add_argument('--max-steps', type=int, default=200,
                        help='Maximum optimization steps (default: 200)')
    parser.add_argument('--pressure', type=float, default=0.0,
                        help='Target pressure in GPa (only for mode 3, default: 0.0)')
    
    args = parser.parse_args()
    
    if args.mode == 3 and args.pressure <= 0:
        parser.error("--pressure must be > 0 when using mode 3 (high pressure)")
    
    os.makedirs(args.relaxed_dir, exist_ok=True)
    optimize_structure(args.cif_file, args.relaxed_dir, 
                      args.mode, args.fmax, args.max_steps, args.pressure)


if __name__ == "__main__":
    main()
