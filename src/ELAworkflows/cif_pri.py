from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import os
import glob

def cif_to_primitive_poscar(cif_path, poscar_path=None):
    structure = Structure.from_file(cif_path)
    
    analyzer = SpacegroupAnalyzer(structure, symprec=0.01)
    primitive_structure = analyzer.get_primitive_standard_structure()
    primitive_structure.to(fmt="poscar", filename=poscar_path)

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Convert CIF to primitive POSCAR')
    parser.add_argument('--input', '-i', help='Input CIF file (default: all *.cif)')
    parser.add_argument('--output', '-o', default='primitive_POSCAR', help='Output POSCAR file')
    args = parser.parse_args()
    
    if args.input:
        cif_files = [args.input]
    else:
        cif_files = glob.glob("*.cif")
    
    if not cif_files:
        print("No CIF files found in current directory.")
    else:
        for cif_file in cif_files:
            cif_to_primitive_poscar(cif_file, args.output)
            print(f"Converted {cif_file} to {args.output}")


if __name__ == "__main__":
    main()
