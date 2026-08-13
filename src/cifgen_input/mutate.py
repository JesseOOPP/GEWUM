from ase.io import read, write
from ase import Atoms
from ase.geometry import get_distances
import os
import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.spatial.distance import pdist, squareform
import shutil

def check_min_distance(structure, min_distance):
    if len(structure) < 2:
        return True
    
    distances = structure.get_all_distances()
    np.fill_diagonal(distances, np.inf) 
    
    min_dist = np.min(distances)
    return min_dist >= min_distance

def apply_lattice_strain(structure, strain_factor):
    new_structure = structure.copy()
    
    strain_tensor = np.random.uniform(-strain_factor, strain_factor, (3, 3))
    strain_tensor = (strain_tensor + strain_tensor.T) / 2 
    
    cell = new_structure.get_cell()
    
    new_cell = cell @ (np.eye(3) + strain_tensor)
    new_structure.set_cell(new_cell, scale_atoms=True)
    
    return new_structure

def rotate_atoms(structure, rotation_probability, max_rotation_angle):
    new_structure = structure.copy()
    
    for i in range(len(new_structure)):
        if np.random.random() < rotation_probability:
            rotation_angles = np.random.uniform(-max_rotation_angle, max_rotation_angle, 3)
            rotation = R.from_euler('xyz', rotation_angles, degrees=True)
            
            pos = new_structure.positions[i]
            
            cell_center = new_structure.get_cell().sum(axis=0) / 2
            relative_pos = pos - cell_center
            
            rotated_pos = rotation.apply(relative_pos)
            
            new_structure.positions[i] = cell_center + rotated_pos
    
    return new_structure

def calculate_structure_fingerprint(structure, tolerance=0.1):
    cell_params = structure.cell.cellpar()
    
    scaled_positions = structure.get_scaled_positions()
    discrete_positions = np.round(scaled_positions / tolerance) * tolerance
    
    features = np.concatenate([
        cell_params,
        discrete_positions.flatten()
    ])
    
    return tuple(np.round(features, 6))

def is_structure_unique(new_structure, existing_fingerprints, tolerance=0.1):
    new_fingerprint = calculate_structure_fingerprint(new_structure, tolerance)
    
    for existing_fp in existing_fingerprints:
        if len(new_fingerprint) != len(existing_fp):
            continue
            
        is_similar = all(abs(a - b) < tolerance for a, b in zip(new_fingerprint, existing_fp))
        if is_similar:
            return False
    
    return True

def generate_modified_structure(original_structure, params, existing_fingerprints, tolerance=0.1):
    max_attempts = 200  
    
    for attempt in range(max_attempts):
        new_structure = original_structure.copy()
        
        if np.random.random() < params['lattice_strain_probability']:
            new_structure = apply_lattice_strain(new_structure, params['max_strain'])
        
        random_stdev = np.random.uniform(params['stdev_min'], params['stdev_max'])
        random_seed = np.random.randint(1, 100000)
        rng = np.random.RandomState(random_seed)
        new_structure.rattle(stdev=random_stdev, rng=rng)
        
        if np.random.random() < params['atom_rotation_probability']:
            new_structure = rotate_atoms(new_structure, 
                                       params['rotation_probability_per_atom'],
                                       params['max_rotation_angle'])
        
        if (check_min_distance(new_structure, params['min_distance']) and
            is_structure_unique(new_structure, existing_fingerprints, tolerance)):
            return new_structure
    
    print(f"Warning: Could not generate unique structure satisfying minimum distance after {max_attempts} attempts")
    return None

def read_params_from_file(filename='INPUT'):
    params = {
        'input_file': 'input.cif', 
        'num_structures': 100,
        'min_distance': 1.0,
        'stdev_min': 0.2,
        'stdev_max': 0.5,
        'lattice_strain_probability': 0.4,
        'max_strain': 0.03,
        'atom_rotation_probability': 0.3,
        'rotation_probability_per_atom': 0.15,
        'max_rotation_angle': 20.0,
        'similarity_tolerance': 0.1  
    }
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('!'):
                    continue
                
                line = line.split('#')[0].split('!')[0].strip()
                
                if '=' in line:
                    key, value = line.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    if key == 'input_file':
                        params[key] = value 
                    else:
                        try:
                            if '.' in value:
                                params[key] = float(value)
                            else:
                                params[key] = int(value)
                        except ValueError:
                            params[key] = value
    
    except FileNotFoundError:
        print(f"Warning: INPUT file '{filename}' not found. Using default parameters.")
    
    return params

def main():
    params = read_params_from_file('INPUT')
    
    print("Parameters used:")
    for key, value in params.items():
        print(f"  {key}: {value}")
    
    input_file = params.get('input_file', 'input.cif')
    
    try:
        original_structure = read(input_file)
        
        distances = original_structure.get_all_distances()
        np.fill_diagonal(distances, np.inf)
        original_min_dist = np.min(distances)
        print(f"Input file: {input_file}")
        print(f"Min_dist: {original_min_dist:.3f}")
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        print("Please make sure the file exists or specify a different file in the INPUT file.")
        return
    except Exception as e:
        print(f"Error reading file '{input_file}': {e}")
        return
    
    output_dir = 'mutated_structures'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    successful_structures = 0
    duplicate_count = 0
    existing_fingerprints = []
    
    original_fingerprint = calculate_structure_fingerprint(original_structure, params.get('similarity_tolerance', 0.1))
    existing_fingerprints.append(original_fingerprint)
    
    i = 0
    while successful_structures < params['num_structures'] and i < params['num_structures'] * 3:  
        i += 1
        
        modified_structure = generate_modified_structure(
            original_structure, 
            params, 
            existing_fingerprints,
            params.get('similarity_tolerance', 0.1)
        )
        
        if modified_structure is not None:
            new_fingerprint = calculate_structure_fingerprint(
                modified_structure, 
                params.get('similarity_tolerance', 0.1)
            )
            existing_fingerprints.append(new_fingerprint)
            
            output_filename = f'{output_dir}/mutated_structure_{successful_structures+1}.cif'
            write(output_filename, modified_structure)
            successful_structures += 1
            
            distances = modified_structure.get_all_distances()
            np.fill_diagonal(distances, np.inf)
            actual_min_dist = np.min(distances)
            print(f"Generated structure {successful_structures}: min_dist = {actual_min_dist:.3f}")
        else:
            duplicate_count += 1
            print(f"Attempt {i}: Duplicate structure or failed to generate")
    
    original_output_filename = f'{output_dir}/original_structure.cif'
    write(original_output_filename, original_structure)
    print(f"Original structure saved as: {original_output_filename}")
    
    print(f"\nAll Done with {successful_structures}/{params['num_structures']} unique files")
    print(f"Duplicates skipped: {duplicate_count}")
    print(f"Total attempts: {i}")
    print(f"Total files in output directory: {successful_structures + 1} (including original structure)")

if __name__ == "__main__":
    main()
