#!/usr/bin/env python
"""
Single CIF file energy calculation module.
Used by calc_energy.sh for parallel processing.
"""
import os
import sys
import csv
from ase.io import read
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from mattersim.forcefield import MatterSimCalculator


def get_chemical_formula(atoms):
    """Get chemical formula from ASE atoms object"""
    symbols = atoms.get_chemical_symbols()
    from collections import Counter
    counts = Counter(symbols)
    
    formula = ""
    for elem in sorted(counts.keys()):
        count = counts[elem]
        if count == 1:
            formula += elem
        else:
            formula += f"{elem}{count}"
    return formula


def get_space_group(cif_file):
    """Get space group from CIF file"""
    try:
        struct = Structure.from_file(cif_file)
        analyzer = SpacegroupAnalyzer(struct)
        return str(analyzer.get_space_group_number())
    except Exception:
        return "Unknown"


def calculate_single_energy(cif_file, output_csv, device="cpu"):
    """
    Calculate energy for a single CIF file and write to CSV.
    
    Args:
        cif_file: Path to CIF file
        output_csv: Path to output CSV file
        device: Device for MatterSim (cpu/cuda)
    """
    if not os.path.exists(cif_file):
        print(f"Error: CIF file not found: {cif_file}", file=sys.stderr)
        return False
    
    try:
        atoms = read(cif_file)
        calc = MatterSimCalculator(device=device)
        atoms.calc = calc
        
        total_energy = atoms.get_potential_energy()
        num_atoms = len(atoms)
        energy_per_atom = total_energy / num_atoms
        
        chemical_formula = get_chemical_formula(atoms)
        cif_base_name = os.path.splitext(os.path.basename(cif_file))[0]
        space_group = get_space_group(cif_file)
        cif_path = os.path.abspath(cif_file)
        
        with open(output_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                chemical_formula,
                cif_base_name,
                f"{total_energy:.6f}",
                f"{energy_per_atom:.6f}",
                cif_path,
                space_group
            ])
        
        print(f"Processed: {cif_file} -> {energy_per_atom:.6f} eV/atom")
        return True
        
    except Exception as e:
        print(f"Error processing {cif_file}: {e}", file=sys.stderr)
        return False


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Single CIF energy calculation')
    parser.add_argument('cif_file', help='Path to CIF file')
    parser.add_argument('output_csv', help='Output CSV file path')
    parser.add_argument('--device', default='cpu', help='Device (cpu/cuda)')
    args = parser.parse_args()
    
    success = calculate_single_energy(args.cif_file, args.output_csv, args.device)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
