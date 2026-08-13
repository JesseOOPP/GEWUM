"""
GEWUM Chiral Structure Classification Module
Classifies CIF structures by their space group symmetry type
"""
import os
from shutil import move
from pymatgen.io.cif import CifParser

CHIRAL_SYMMETRY_NUMBERS = [
    3, 4, 5, 16, 17, 18, 19, 20, 21, 22, 23, 24, 75, 76, 77, 78, 79, 80, 
    89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 143, 144, 145, 146, 149, 150, 
    151, 152, 153, 154, 155, 168, 169, 170, 171, 172, 173, 177, 178, 179, 
    180, 181, 182, 195, 196, 197, 198, 199, 207, 208, 209, 210, 211, 212, 
    213, 214
]

P1_SYMMETRY_NUMBERS = [1]

NONCENTRO_SYMMETRY_NUMBERS = [
    6, 7, 8, 9, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 
    39, 40, 41, 42, 43, 44, 45, 46, 81, 82, 99, 100, 101, 102, 103, 104, 
    105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 
    119, 120, 121, 122, 156, 157, 158, 159, 160, 161, 174, 183, 184, 185, 
    186, 187, 188, 189, 190, 215, 216, 217, 218, 219, 220
]

INVERSION_SYMMETRY_NUMBERS = [
    2, 10, 11, 12, 13, 14, 15, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 
    58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 83, 
    84, 85, 86, 87, 88, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 
    133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 147, 148, 162, 163, 
    164, 165, 166, 167, 175, 176, 191, 192, 193, 194, 200, 201, 202, 203, 
    204, 205, 206, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230
]


def classify_space_group(space_group_number):
    """
    Classify a space group number into symmetry category.
    
    Args:
        space_group_number: International space group number (1-230)
    
    Returns:
        Category name: 'chiral', 'P1', 'noncentro', 'inversion', or None
    """
    if space_group_number in CHIRAL_SYMMETRY_NUMBERS:
        return 'chiral'
    elif space_group_number in P1_SYMMETRY_NUMBERS:
        return 'P1'
    elif space_group_number in NONCENTRO_SYMMETRY_NUMBERS:
        return 'noncentro'
    elif space_group_number in INVERSION_SYMMETRY_NUMBERS:
        return 'inversion'
    return None


def main(work_dir='.'):
    """
    Main function to classify and sort CIF files by symmetry type.
    
    Args:
        work_dir: Working directory containing CIF files
    """
    import argparse
    parser = argparse.ArgumentParser(description='GEWUM Chiral Structure Classification')
    parser.add_argument('--dir', '-d', default='.',
                        help='Working directory containing CIF files (default: current directory)')
    args = parser.parse_args()
    work_dir = args.dir
    
    categories = ['chiral', 'P1', 'noncentro', 'inversion']
    for category in categories:
        os.makedirs(os.path.join(work_dir, category), exist_ok=True)
    
    cif_files = [f for f in os.listdir(work_dir) if f.endswith('.cif')]
    
    log_path = os.path.join(work_dir, 'classification_log.txt')
    with open(log_path, 'w') as log_file:
        for cif_file in cif_files:
            cif_path = os.path.join(work_dir, cif_file)
            try:
                parser = CifParser(cif_path)
                structure = parser.parse_structures(primitive=True)[0]
                
                space_group_number = structure.get_space_group_info()[1]
                target_dir = classify_space_group(space_group_number)
                
                if target_dir:
                    dest_path = os.path.join(work_dir, target_dir, cif_file)
                    move(cif_path, dest_path)
                    log_message = f"Moved {cif_file} to {target_dir} (SG: {space_group_number})"
                else:
                    log_message = f"Kept {cif_file} in current directory (SG: {space_group_number})"
                
                log_file.write(log_message + '\n')
                print(log_message)
            
            except Exception as e:
                error_message = f"Error processing {cif_file}: {e}"
                log_file.write(error_message + '\n')
                print(error_message)
    
    print(f"\nClassification complete. Log saved to: {log_path}")


if __name__ == "__main__":
    main()
