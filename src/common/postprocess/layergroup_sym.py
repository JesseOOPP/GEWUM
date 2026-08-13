"""
GEWUM Layer Group Symmetry Analysis Module
Classifies 2D structures based on layer group symmetry and moves to appropriate directories
"""
import os
import re
from shutil import move
import numpy as np
from ase.io import read
import spglib


# Non-centrosymmetric layer group numbers for 2D materials
NONCENTRO_LAYER_GROUPS = [
    3, 4, 5, 8, 9, 10, 11, 12, 14, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 
    29, 30, 31, 32, 33, 34, 35, 36, 49, 50, 53, 54, 55, 56, 57, 58, 59, 60, 
    65, 67, 68, 69, 70, 73, 74, 76, 77, 78, 79
]


def get_simplified_formula(cif_file):
    """
    Extract and simplify chemical formula from CIF file.
    
    Args:
        cif_file: Path to CIF file
    
    Returns:
        Simplified chemical formula string
    """
    with open(cif_file, 'r') as file:
        content = file.read()
    
    formula_match = re.search(r'_chemical_formula_structural\s+(.*)', content)
    if not formula_match:
        formula_match = re.search(r'_chemical_formula_sum\s+(.*)', content)
    
    if not formula_match:
        return os.path.splitext(os.path.basename(cif_file))[0]
    
    chemical_formula = formula_match.group(1).strip().strip("'\"")
    chemical_formula = chemical_formula.replace(' ', '').replace('.', '')
    
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', chemical_formula)
    counts = []
    for elem, count_str in elements:
        count = int(count_str) if count_str else 1
        counts.append(count)
    
    gcd = np.gcd.reduce(counts) if len(counts) > 0 else 1
    
    simplified_formula = ""
    for elem, count_str in elements:
        count = int(count_str) if count_str else 1
        new_count = count // gcd
        if new_count > 1:
            simplified_formula += f"{elem}{new_count}"
        else:
            simplified_formula += elem
    
    illegal_chars = r'[\\/:*?"<>|\s]'
    simplified_formula = re.sub(illegal_chars, '_', simplified_formula)
    simplified_formula = re.sub(r'_+', '_', simplified_formula).strip('_')
    
    return simplified_formula


def get_layer_group_number(cif_file, symprec=0.05):
    """
    Get layer group number for a 2D structure.
    
    Args:
        cif_file: Path to CIF file
        symprec: Symmetry precision for analysis
    
    Returns:
        Layer group number
    """
    ase_structure = read(cif_file)
    
    lattice = ase_structure.get_cell()
    positions = ase_structure.get_scaled_positions()
    numbers = ase_structure.get_atomic_numbers()
    
    cell = (lattice, positions, numbers)
    
    layergroup_result = spglib.get_layergroup(cell, symprec=symprec)
    return layergroup_result.number


def classify_and_move_cif(cif_file, work_dir, symprec=0.05):
    """
    Classify a 2D CIF file based on layer group symmetry and move to appropriate directory.
    
    Args:
        cif_file: Path to CIF file
        work_dir: Working directory containing output subdirectories
        symprec: Symmetry precision for analysis
    
    Returns:
        Tuple of (success, message)
    """
    try:
        layergroup_number = get_layer_group_number(cif_file, symprec)
        simplified_formula = get_simplified_formula(cif_file)
        
        if layergroup_number in NONCENTRO_LAYER_GROUPS:
            target_dir = os.path.join(work_dir, 'noncentro')
        else:
            target_dir = os.path.join(work_dir, 'others')
        
        os.makedirs(target_dir, exist_ok=True)
        
        base_filename = f"{layergroup_number}_{simplified_formula}"
        new_filename = f"{base_filename}.cif"
        
        counter = 1
        while os.path.exists(os.path.join(target_dir, new_filename)):
            new_filename = f"{base_filename}_{counter}.cif"
            counter += 1
        
        move(cif_file, os.path.join(target_dir, new_filename))
        
        message = f"Renamed and moved {os.path.basename(cif_file)} to {os.path.basename(target_dir)}/{new_filename} (Layer group: {layergroup_number})"
        return True, message
        
    except Exception as e:
        message = f"Error processing {os.path.basename(cif_file)}: {e}"
        return False, message


def main(work_dir=None, symprec=0.05):
    """
    Main function to classify 2D CIF files based on layer group symmetry.
    
    Args:
        work_dir: Working directory containing CIF files (defaults to current directory)
        symprec: Symmetry precision for analysis
    """
    import argparse
    parser = argparse.ArgumentParser(description='GEWUM Layer Group Symmetry Analysis')
    parser.add_argument('--dir', '-d', default=None,
                        help='Working directory containing CIF files (default: current directory)')
    parser.add_argument('--symprec', '-s', type=float, default=0.05,
                        help='Symmetry precision for analysis (default: 0.05)')
    args = parser.parse_args()
    
    if args.dir:
        work_dir = args.dir
    elif work_dir is None:
        work_dir = os.getcwd()
    symprec = args.symprec
    
    noncentro_dir = os.path.join(work_dir, 'noncentro')
    others_dir = os.path.join(work_dir, 'others')
    os.makedirs(noncentro_dir, exist_ok=True)
    os.makedirs(others_dir, exist_ok=True)
    
    cif_files = [f for f in os.listdir(work_dir) if f.endswith('.cif')]
    
    if not cif_files:
        print("No CIF files found in the working directory.")
        return
    
    log_path = os.path.join(work_dir, 'classification_log.txt')
    success_count = 0
    error_count = 0
    
    with open(log_path, 'w') as log_file:
        for cif_file in cif_files:
            cif_path = os.path.join(work_dir, cif_file)
            success, message = classify_and_move_cif(cif_path, work_dir, symprec)
            
            log_file.write(message + '\n')
            print(message)
            
            if success:
                success_count += 1
            else:
                error_count += 1
    
    print(f"\nProcessing summary:")
    print(f"Successfully processed: {success_count} files")
    print(f"Files with errors: {error_count} files")
    print(f"Log saved to: {log_path}")


if __name__ == "__main__":
    main()
