import os
import yaml
import re
import numpy as np
from ase.io import read, write
from itertools import product

ELEMENT_GROUPS = {
    "G0": ["H"],
    "G1": ["Li", "Na", "K", "Rb", "Cs"],
    "G2": ["Be", "Mg", "Ca", "Sr", "Ba"],
    "G3": ["B", "Al", "Ga", "In", "Tl"],
    "G4": ["C", "Si", "Ge", "Sn", "Pb"],
    "G5": ["N", "P", "As", "Sb","Bi"],
    "G6-1": ["O"],
    "G6-2": ["S", "Se", "Te"],
    "G7": ["F", "Cl", "Br", "I"],
    "G8": ["He", "Ne", "Ar", "Kr", "Xe"],
    "T1": ["Sc", "Y", "Ti", "Zr", "Hf", "V", "Nb", "Ta", "Cr", "Mo", "W", "Mn", "Re", "Fe", "Ru", "Os", "Co", "Rh", "Ir", "Ni", "Pd", "Pt"],
    "T2": ["Cu", "Ag", "Au", "Zn", "Cd", "Hg"],
    "LAN": ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"],
    "ACT": ["Ac", "Th", "Pa", "U", "Np", "Pu"]
}

def clean_yaml_content(content):
    """Remove null bytes and other problematic characters from YAML content"""
    cleaned = content.replace('\x00', '')
    cleaned = re.sub(r'[^\x20-\x7E\t\n\r]', '', cleaned)
    return cleaned

def expand_group_references(replacements):
    expanded = {}
    for element, replacement in replacements.items():
        if isinstance(replacement, str) and replacement in ELEMENT_GROUPS:
            expanded[element] = ELEMENT_GROUPS[replacement]
        elif isinstance(replacement, list):
            expanded[element] = replacement
        else:
            raise ValueError(f"Invalid replacement definition: {element} -> {replacement}")
    return expanded

def generate_combinations(replacements):
    keys = list(replacements.keys())
    value_lists = [replacements[key] for key in keys]
    for combination in product(*value_lists):
        yield dict(zip(keys, combination))

def process_cif(cif_path, replacements, output_dir):
    if not os.path.exists(cif_path):
        print(f"Warning: CIF file {cif_path} not found. Skipping.")
        return
    
    try:
        structure = read(cif_path)
    except Exception as e:
        print(f"Error reading {cif_path}: {str(e)}")
        return

    # Validate cell integrity to avoid "scaling factors" errors downstream
    cell = np.array(structure.get_cell())
    if cell.ndim != 2 or cell.shape != (3, 3) or np.linalg.matrix_rank(cell) < 3:
        print(f"Error: Invalid or degenerate cell in {cif_path}. Skipping.")
        return

    # Cache cell and scaled positions - needed to re-anchor after copy
    original_cell = cell.copy()
    original_scaled = structure.get_scaled_positions()

    base_name = os.path.splitext(os.path.basename(cif_path))[0]
    
    try:
        expanded_replacements = expand_group_references(replacements)
    except ValueError as e:
        print(f"Error in replacement rules for {cif_path}: {str(e)}")
        return
    
    present_elements = set(atom.symbol for atom in structure)
    valid_replacements = {k: v for k, v in expanded_replacements.items() if k in present_elements}
    
    if not valid_replacements:
        print(f"Warning: No replaceable elements found in {cif_path}. Skipping.")
        return
    
    for i, combo in enumerate(generate_combinations(valid_replacements)):
        new_structure = structure.copy()
        for idx, atom in enumerate(new_structure):
            if atom.symbol in combo:
                new_structure[idx].symbol = combo[atom.symbol]

        # Explicitly re-anchor cell & positions to avoid ASE internal
        # "The number of scaling factors must be 1 or 3" during CIF write
        new_structure.set_cell(original_cell, scale_atoms=False)
        new_structure.set_scaled_positions(original_scaled)
        
        output_path = os.path.join(output_dir, f"{base_name}_combo_{i+1}.cif")
        try:
            write(output_path, new_structure, format='cif')
            print(f"Generated {output_path} with replacements: {combo}")
        except Exception as e:
            print(f"Error writing {output_path}: {str(e)}")

def main():
    input_file = "replacements.yaml"
    output_dir = "substituted_structures"
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        with open(input_file, 'r') as f:
            raw_content = f.read()
            cleaned_content = clean_yaml_content(raw_content)
            config = yaml.safe_load(cleaned_content) or {}
            
            if not isinstance(config, dict):
                raise ValueError("Invalid YAML structure. Expected a dictionary of CIF files and their replacements.")
                
            print(f"Successfully loaded configuration for {len(config)} CIF files")
            
    except FileNotFoundError:
        print(f"Error: Configuration file {input_file} not found")
        return
    except yaml.YAMLError as e:
        print(f"YAML parsing error: {str(e)}")
        print("Problematic content:")
        print(cleaned_content[:200] + "..." if len(cleaned_content) > 200 else cleaned_content)
        return
    except Exception as e:
        print(f"Unexpected error: {str(e)}")
        return
    
    for cif_file, replacements in config.items():
        print(f"\nProcessing {cif_file} with replacements: {replacements}")
        process_cif(cif_file, replacements, output_dir)

if __name__ == "__main__":
    main()
