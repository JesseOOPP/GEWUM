"""
Supercell Generation Script for PT Workflow
Generate supercells from input CIF files using pymatgen
"""
from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter
import os
import sys
import argparse


def make_supercell(input_file, scaling_matrix, output_file):
    """
    Generate supercell from input CIF file
    
    Args:
        input_file: Path to input CIF file
        scaling_matrix: List of 3 integers or 3x3 matrix for supercell scaling
        output_file: Path to output CIF file
    """
    try:
        structure = Structure.from_file(input_file)
        original_num_atoms = len(structure)
        
        structure.make_supercell(scaling_matrix)
        
        writer = CifWriter(structure)
        writer.write_file(output_file)
        
        print(f"[OK] Supercell generated: {output_file}")
        print(f"  Original: {original_num_atoms} atoms")
        print(f"  Supercell: {len(structure)} atoms")
        print(f"  Scaling: {scaling_matrix}")
        
        return True
        
    except Exception as e:
        print(f"[FAIL] Error generating supercell from {input_file}: {e}")
        return False


def batch_supercell(input_dir='.', scaling_matrix=[2, 2, 1], output_dir='supercell_structures'):
    """
    Batch generate supercells for all CIF files in a directory
    
    Args:
        input_dir: Directory containing input CIF files
        scaling_matrix: Scaling matrix for supercell
        output_dir: Output directory for supercell files
    """
    os.makedirs(output_dir, exist_ok=True)
    
    cif_files = [f for f in os.listdir(input_dir) if f.endswith('.cif')]
    
    if not cif_files:
        print(f"No CIF files found in {input_dir}")
        return
    
    print(f"Found {len(cif_files)} CIF file(s)")
    print(f"Scaling matrix: {scaling_matrix}")
    print(f"Output directory: {output_dir}\n")
    
    success_count = 0
    for cif_file in cif_files:
        input_path = os.path.join(input_dir, cif_file)
        output_filename = f"supercell_{cif_file}"
        output_path = os.path.join(output_dir, output_filename)
        
        if make_supercell(input_path, scaling_matrix, output_path):
            success_count += 1
    
    print(f"\n{'='*60}")
    print(f"Supercell generation complete: {success_count}/{len(cif_files)} files processed")
    print(f"{'='*60}")


def main():
    """Main function with command line arguments"""
    parser = argparse.ArgumentParser(
        description='Generate supercells from CIF files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file (default: 2x2x2)
  python supercell.py -i input.cif
  
  # Single file with custom matrix
  python supercell.py -i input.cif --matrix 2 2 1 -o supercell.cif
  
  # Batch mode (process all CIF files)
  python supercell.py --batch --matrix 3 3 1
        """
    )
    
    parser.add_argument('-i', '--input', help='Input CIF file')
    parser.add_argument('-o', '--output', help='Output CIF file (default: supercell_<input>.cif)')
    parser.add_argument('--matrix', nargs=3, type=int, default=[2, 2, 2], help='Scaling matrix (3 integers, default: 2 2 2)')
    parser.add_argument('--batch', action='store_true', help='Batch process all CIF files in current directory')
    parser.add_argument('--output-dir', default='supercell_structures', help='Output directory for batch mode (default: supercell_structures)')
    
    args = parser.parse_args()
    
    if args.batch:
        print("=" * 60)
        print("Supercell Generation (Batch Mode)")
        print("=" * 60)
        batch_supercell('.', args.matrix, args.output_dir)
        
    else:
        if not args.input:
            print("Error: -i/--input is required (or use --batch for batch mode)")
            parser.print_help()
            sys.exit(1)
        
        if not args.output:
            base_name = os.path.splitext(args.input)[0]
            args.output = f"supercell_{base_name}.cif"
        
        print("=" * 60)
        print("Supercell Generation (Single File Mode)")
        print("=" * 60)
        
        make_supercell(args.input, args.matrix, args.output)


if __name__ == '__main__':
    main()
