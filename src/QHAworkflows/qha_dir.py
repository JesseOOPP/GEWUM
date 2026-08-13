import os
import shutil

def get_file_prefix(filename):
    if filename.endswith('.cif'):
        return filename[:-4] 
    return None

def create_and_copy_cif_files():
    cif_files = [f for f in os.listdir('.') if f.endswith('.cif')]
    
    for cif_file in cif_files:
        prefix = get_file_prefix(cif_file)

        folder_name = prefix
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        
        src_path = os.path.join('.', cif_file)
        dest_path = os.path.join('.', folder_name, cif_file)
        shutil.copy(src_path, dest_path)

if __name__ == "__main__":
    create_and_copy_cif_files()
    print("All Done")
