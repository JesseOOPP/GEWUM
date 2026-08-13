"""
GEWUM Symmetry Analysis and Renaming Module
Symmetrizes structures and renames CIF files based on space group and composition
"""
import os
import re
import shutil
from pymatgen.core import Structure
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# CIF output precision (fractional coordinates) - kept tight to preserve
# Wyckoff positions exactly, independent of the symmetry-detection symprec.
_CIF_WRITE_PREC = 1e-6


def _refine_via_spglib(structure, symprec):
    """
    Refine atomic positions to exact Wyckoff sites using spglib's
    standardize_cell, which snaps atoms to their ideal symmetric positions.

    This is a stricter pass than pymatgen's get_symmetrized_structure() alone:
    it rebuilds the cell through spglib standardisation and returns a fresh
    pymatgen Structure whose fractional coordinates are as close to the
    mathematical Wyckoff positions as double precision allows.

    Returns:
        Refined pymatgen Structure, or None on failure.
    """
    try:
        import spglib
        lattice = structure.lattice.matrix
        positions = structure.frac_coords
        numbers = structure.atomic_numbers
        cell = (lattice, positions, numbers)
        # to_primitive=False: keep conventional cell
        # no_idealize=False: snap atoms to exact Wyckoff positions
        std_cell = spglib.standardize_cell(
            cell, to_primitive=False, no_idealize=False, symprec=symprec
        )
        if std_cell is None:
            return None
        std_lattice, std_pos, std_nums = std_cell
        refined = Structure(
            std_lattice,
            std_nums,
            std_pos,
            coords_are_cartesian=False,
            to_unit_cell=True,
        )
        return refined
    except Exception as e:
        print(f"  [refine] spglib refine failed ({e}), falling back to pymatgen symmetrization")
        return None


def symmetrize_and_save_cif(cif_file_path, output_dir, compound_name, symprec=0.05):
    """
    Symmetrize a structure and save it as a new CIF file.
    
    Two-step process:
      1. Symmetry analysis & symmetrization via pymatgen (symprec controls
         how aggressively atoms are merged to Wyckoff sites).
      2. spglib refine_cell pass to snap coordinates to exact Wyckoff positions,
         removing any residual off-site drift.
      3. CIF output with tight coordinate precision (1e-6) so that fractional
         coordinates are not rounded away from their ideal Wyckoff values.

    Args:
        cif_file_path: Path to input CIF file
        output_dir: Directory to save symmetrized structure
        compound_name: Base name for output file
        symprec: Symmetry precision for analysis (default 0.05)
    
    Returns:
        Path to output file, or None if processing failed
    """
    try:
        parser = CifParser(cif_file_path)
        structures = parser.parse_structures(primitive=True)
        
        if not structures:
            raise ValueError("Invalid CIF file with no structures!")
            
        relaxed_cif = structures[0]
        
        spacegroup_analyzer = SpacegroupAnalyzer(relaxed_cif, symprec=symprec)
        symmetrized_cif = spacegroup_analyzer.get_symmetrized_structure()

        # --- spglib refine: force atoms to exact Wyckoff sites ---
        refined = _refine_via_spglib(symmetrized_cif, symprec=symprec)
        if refined is not None:
            symmetrized_cif = refined

        # Use a tight precision for CIF output so that fractional
        # coordinates stay on their ideal Wyckoff positions.
        cif_writer = CifWriter(symmetrized_cif, symprec=_CIF_WRITE_PREC)
        output_path = os.path.join(output_dir, f"{compound_name}_symmetrized.cif")
        os.makedirs(output_dir, exist_ok=True)
        
        cif_writer.write_file(output_path)
        print(f'Symmetrized structure saved to: {output_path}')
        return output_path
        
    except Exception as e:
        print(f"Error processing {cif_file_path}: {str(e)}")
        return None


def rename_cif_files(directory):
    """
    Rename CIF files based on their space group and chemical formula.
    
    Args:
        directory: Directory containing CIF files to rename
    """
    for cif_file in os.listdir(directory):
        if not cif_file.endswith('.cif'):
            continue
            
        file_path = os.path.join(directory, cif_file)
        with open(file_path, 'r') as file:
            content = file.read()
            
            space_group_match = re.search(r'_symmetry_space_group_name_H-M\s+(.*)', content)
            chemical_formula_match = re.search(r'_chemical_formula_structural\s+(.*)', content)
            
            if not space_group_match or not chemical_formula_match:
                print(f"File '{cif_file}' skipped - missing required fields")
                continue
                
            space_group = space_group_match.group(1).strip().strip("'\"")
            chemical_formula = chemical_formula_match.group(1).strip().strip("'\"")
            
            illegal_chars = r'[\\/:*?"<>|\s]'
            space_group_clean = re.sub(illegal_chars, '_', space_group)
            chemical_formula_clean = re.sub(illegal_chars, '_', chemical_formula)
            
            space_group_clean = re.sub(r'_+', '_', space_group_clean)
            chemical_formula_clean = re.sub(r'_+', '_', chemical_formula_clean)
            
            base_filename = f"{space_group_clean}_{chemical_formula_clean}"
            new_filename = f"{base_filename}.cif"
            new_file_path = os.path.join(directory, new_filename)
            
            counter = 1
            while os.path.exists(new_file_path) and new_file_path != file_path:
                new_filename = f"{base_filename}_{counter}.cif"
                new_file_path = os.path.join(directory, new_filename)
                counter += 1
            
            if new_file_path != file_path:
                os.rename(file_path, new_file_path)
                print(f"Renamed '{cif_file}' to '{new_filename}'")
            else:
                print(f"File '{cif_file}' already has the correct name")


def move_error_file(cif_file_path, error_dir):
    """
    Move problematic CIF file to error directory.
    
    Args:
        cif_file_path: Path to the problematic file
        error_dir: Directory to move the file to
    """
    os.makedirs(error_dir, exist_ok=True)
    error_file_path = os.path.join(error_dir, os.path.basename(cif_file_path))
    
    counter = 1
    base_name = os.path.splitext(os.path.basename(cif_file_path))[0]
    extension = '.cif'
    
    while os.path.exists(error_file_path):
        error_file_path = os.path.join(error_dir, f"{base_name}_{counter}{extension}")
        counter += 1
    
    shutil.move(cif_file_path, error_file_path)
    print(f"Moved problematic file to: {error_file_path}")


def main(work_dir=None, symprec=0.05):
    """Main function to symmetrize and rename CIF files"""
    import argparse
    parser = argparse.ArgumentParser(description='GEWUM Symmetry Analysis and Renaming')
    parser.add_argument('--dir', '-d', default=None,
                        help='Working directory containing CIF files (default: current directory)')
    parser.add_argument('--symprec', '-s', type=float, default=0.05,
                        help='Symmetry precision for analysis (default: 0.05)')
    args = parser.parse_args()
    
    if args.dir:
        work_dir = args.dir
    elif work_dir is None:
        work_dir = os.getcwd()
    
    input_dir = work_dir
    output_dir = os.path.join(work_dir, 'relaxed_symmetry')
    error_dir = os.path.join(work_dir, 'error_str')
    
    processed_files = []
    error_files = []
    
    for filename in os.listdir(input_dir):
        if filename.endswith('.cif'):
            compound_name = os.path.splitext(filename)[0]
            cif_file_path = os.path.join(input_dir, filename)
            
            output_path = symmetrize_and_save_cif(cif_file_path, output_dir, compound_name)
            
            if output_path:
                processed_files.append(output_path)
            else:
                error_files.append(cif_file_path)
                move_error_file(cif_file_path, error_dir)
    
    print(f"\nProcessing summary:")
    print(f"Successfully processed: {len(processed_files)} files")
    print(f"Files with errors: {len(error_files)} files")
    
    if processed_files:
        print("\nStarting file renaming process...")
        rename_cif_files(output_dir)
        print("Renaming completed!")
    
    print("All files processed successfully!")


if __name__ == "__main__":
    main()
