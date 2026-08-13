"""
GEWUM Doping Module
Generate doped structures by randomly replacing specified number of atoms
Uses symmetry analysis to eliminate equivalent site combinations
"""
import os
import sys
import yaml
import re
import random
from itertools import combinations, product
from ase.io import read, write
import numpy as np

try:
    from spglib import get_symmetry_dataset
    SPGLIB_AVAILABLE = True
except ImportError:
    SPGLIB_AVAILABLE = False
    print("Warning: spglib not available. Symmetry optimization disabled.")


def clean_yaml_content(content):
    """Remove null bytes and other problematic characters from YAML content"""
    cleaned = content.replace('\x00', '')
    cleaned = re.sub(r'[^\x20-\x7E\t\n\r]', '', cleaned)
    return cleaned


def parse_doping_rule(rule_str):
    """
    Parse doping rule string like 'Li -> Na: 2'
    Returns: (source_element, target_element)
    """
    match = re.match(r'(\w+)\s*->\s*(\w+)', rule_str)
    if not match:
        raise ValueError(f"Invalid doping rule format: {rule_str}")
    return match.group(1), match.group(2)


def get_element_indices(structure, element):
    """Get all indices of a specific element in structure"""
    return [i for i, atom in enumerate(structure) if atom.symbol == element]


def are_configurations_equivalent(structure, indices1, indices2, symprec=1e-5):
    """
    Check if two doping configurations are symmetrically equivalent
    
    Args:
        structure: ASE Atoms object
        indices1, indices2: Two sets of atom indices to be doped
        symprec: Symmetry precision
    
    Returns:
        True if configurations are equivalent under symmetry operations
    """
    if not SPGLIB_AVAILABLE:
        return False
    
    try:
        cell = (structure.get_cell(), 
                structure.get_scaled_positions(), 
                structure.get_atomic_numbers())
        
        symmetry_data = get_symmetry_dataset(cell, symprec=symprec)
        if symmetry_data is None:
            return False
        
        rotations = symmetry_data['rotations']
        translations = symmetry_data['translations']
        
        scaled_pos = structure.get_scaled_positions()
        set1 = set(indices1)
        set2 = set(indices2)
        
        for rot, trans in zip(rotations, translations):
            mapped_set = set()
            
            for idx1 in set1:
                pos1 = scaled_pos[idx1]
                transformed_pos = (np.dot(pos1, rot.T) + trans) % 1.0
                
                for idx2, pos2 in enumerate(scaled_pos):
                    diff = (transformed_pos - pos2) % 1.0
                    diff = np.minimum(diff, 1.0 - diff)
                    if np.all(diff < symprec):
                        mapped_set.add(idx2)
                        break
            
            if mapped_set == set2:
                return True
        
        return False
        
    except Exception as e:
        return False


def filter_symmetry_equivalent_combinations(structure, combinations_list):
    """
    Filter out symmetry-equivalent combinations, keeping only unique ones
    
    Args:
        structure: ASE Atoms object
        combinations_list: List of index tuples
    
    Returns:
        List of symmetry-unique combinations
    """
    if not SPGLIB_AVAILABLE or len(combinations_list) <= 1:
        return combinations_list
    
    unique = []
    for combo in combinations_list:
        is_unique = True
        for existing in unique:
            if are_configurations_equivalent(structure, combo, existing):
                is_unique = False
                break
        if is_unique:
            unique.append(combo)
    
    return unique


def generate_doping_configurations(structure, doping_rules, max_configs=None):
    """
    Generate doped configurations based on rules with symmetry optimization
    
    Args:
        structure: ASE Atoms object
        doping_rules: dict of {source_element: (target_element, count)}
        max_configs: maximum number of configurations to generate (None = all)
    
    Returns:
        List of (doped_structure, description) tuples
    """
    print(f"  Analyzing structure and generating combinations...")
    
    element_all_combos = {}
    element_counts = {}
    
    for source_elem, (target_elem, count) in doping_rules.items():
        indices = get_element_indices(structure, source_elem)
        element_counts[source_elem] = len(indices)
        
        if len(indices) < count:
            print(f"  Warning: Only {len(indices)} {source_elem} atoms found, but {count} requested")
            return []
        
        all_combos = list(combinations(indices, count))
        
        print(f"  {source_elem} ({len(indices)} atoms) -> {target_elem}: replace {count} atoms")
        print(f"    Total possible combinations: {len(all_combos)}")
        
        if SPGLIB_AVAILABLE:
            unique_combos = filter_symmetry_equivalent_combinations(structure, all_combos)
            print(f"    Symmetry-unique combinations: {len(unique_combos)}")
            element_all_combos[source_elem] = unique_combos
        else:
            print(f"    (Symmetry filtering disabled)")
            element_all_combos[source_elem] = all_combos
    
    total_combinations = 1
    for combos in element_all_combos.values():
        total_combinations *= len(combos)
    
    print(f"\n  Total unique configurations: {total_combinations}")
    
    if max_configs is None or max_configs >= total_combinations:
        generate_all = True
        num_to_generate = total_combinations
        print(f"  Generating all {total_combinations} unique configurations")
    else:
        generate_all = False
        num_to_generate = max_configs
        print(f"  Generating {num_to_generate} configurations (sampled from {total_combinations})")
    
    configurations = []
    combo_keys = list(element_all_combos.keys())
    
    if generate_all:
        all_combos = list(product(*[element_all_combos[k] for k in combo_keys]))
    else:
        all_combos = list(product(*[element_all_combos[k] for k in combo_keys]))
        if num_to_generate < len(all_combos):
            all_combos = random.sample(all_combos, num_to_generate)
    
    # Cache original cell and positions to re-anchor after copy
    # (avoids ASE CIF writer "scaling factors" error)
    original_cell = np.array(structure.get_cell())
    original_scaled = structure.get_scaled_positions()
    
    for combo_tuple in all_combos:
        replacement_map = {}
        for i, source_elem in enumerate(combo_keys):
            target_elem = doping_rules[source_elem][0]
            for idx in combo_tuple[i]:
                replacement_map[idx] = target_elem
        
        new_structure = structure.copy()
        for idx, new_elem in replacement_map.items():
            new_structure[idx].symbol = new_elem
        new_structure.set_cell(original_cell, scale_atoms=False)
        new_structure.set_scaled_positions(original_scaled)
        
        desc_parts = []
        for source_elem in combo_keys:
            target_elem = doping_rules[source_elem][0]
            count = doping_rules[source_elem][1]
            desc_parts.append(f"{source_elem}{count}{target_elem}")
        description = "_".join(desc_parts)
        
        configurations.append((new_structure, description))
    
    return configurations


def process_cif_doping(cif_path, doping_config, output_dir):
    """
    Process a single CIF file with doping rules
    
    Args:
        cif_path: Path to input CIF file
        doping_config: dict with doping rules
        output_dir: Output directory
    """
    if not os.path.exists(cif_path):
        print(f"Warning: CIF file {cif_path} not found. Skipping.")
        return
    
    try:
        structure = read(cif_path)
    except Exception as e:
        print(f"Error reading {cif_path}: {str(e)}")
        return
    
    # Validate cell integrity to avoid ASE CIF writer errors
    cell = np.array(structure.get_cell())
    if cell.ndim != 2 or cell.shape != (3, 3) or np.linalg.matrix_rank(cell) < 3:
        print(f"Error: Invalid or degenerate cell in {cif_path}. Skipping.")
        return
    
    base_name = os.path.splitext(os.path.basename(cif_path))[0]
    
    doping_rules = {}
    for rule_key, count in doping_config.items():
        try:
            source_elem, target_elem = parse_doping_rule(rule_key)
            if source_elem in doping_rules:
                prev_target = doping_rules[source_elem][0]
                print(f"Warning: Duplicate doping rule for '{source_elem}' "
                      f"(was '{source_elem} -> {prev_target}', now '{source_elem} -> {target_elem}'). "
                      f"Only the last rule will be used.")
            doping_rules[source_elem] = (target_elem, count)
        except ValueError as e:
            print(f"Error parsing doping rule '{rule_key}': {str(e)}")
            return
    
    print(f"\nProcessing {cif_path}")
    print(f"Doping rules:")
    for rule_key, count in doping_config.items():
        print(f"  {rule_key}: {count}")
    
    configurations = generate_doping_configurations(structure, doping_rules, max_configs=None)
    
    if not configurations:
        print(f"No valid configurations generated for {cif_path}")
        return
    
    total_possible = len(configurations)
    
    # Non-interactive environments (SLURM / pipe): generate all automatically
    if not sys.stdin.isatty():
        print(f"\n  Non-interactive mode: generating all {total_possible} configurations")
    else:
        print(f"\n  Generate all {total_possible} configurations? (y/n): ", end='', flush=True)
        response = input().strip().lower()
        
        if response != 'y':
            print(f"  Enter maximum number of configurations to generate: ", end='', flush=True)
            try:
                max_num = int(input().strip())
                if max_num <= 0:
                    print("  Invalid number. Skipping this file.")
                    return
                if max_num < total_possible:
                    configurations = random.sample(configurations, max_num)
            except (ValueError, EOFError):
                print("  Invalid input. Skipping this file.")
                return
    
    print(f"\n  Generating {len(configurations)} structures...")
    for i, (doped_structure, description) in enumerate(configurations, 1):
        output_filename = f"{base_name}_{description}_{i:03d}.cif"
        output_path = os.path.join(output_dir, output_filename)
        
        try:
            write(output_path, doped_structure, format='cif')
        except Exception as e:
            print(f"  Error writing {output_path}: {str(e)}")
    
    print(f"  [OK] Generated {len(configurations)} doped structures")


def main():
    """Main function for doping workflow"""
    input_file = "doping.yaml"
    output_dir = "doped_structures"
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            raw_content = f.read()
            cleaned_content = clean_yaml_content(raw_content)
            config = yaml.safe_load(cleaned_content) or {}
            
            if not isinstance(config, dict):
                raise ValueError("Invalid YAML structure. Expected a dictionary of CIF files and doping rules.")
            
            print(f"Loaded doping configuration for {len(config)} CIF file(s)")
            
    except FileNotFoundError:
        print(f"Error: Configuration file {input_file} not found")
        print(f"Please run 'gewum PT --mode dp' to copy the template file")
        return
    except yaml.YAMLError as e:
        print(f"YAML parsing error: {str(e)}")
        return
    except Exception as e:
        print(f"Unexpected error loading configuration: {str(e)}")
        return
    
    for cif_file, doping_rules in config.items():
        if not isinstance(doping_rules, dict):
            print(f"Warning: Invalid doping rules for {cif_file}. Skipping.")
            continue
        
        process_cif_doping(cif_file, doping_rules, output_dir)
    
    print(f"\n{'='*60}")
    print(f"Doping workflow completed")
    print(f"Output directory: {output_dir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
