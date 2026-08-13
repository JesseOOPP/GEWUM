import os
import csv
from pymatgen.core import Structure

def main():
    current_dir = os.getcwd()  
    chem_dir_name = os.path.basename(current_dir)
    
    csv_file = os.path.join(current_dir, 'energy_results.csv')
    
    if not os.path.exists(csv_file):
        print(f"Error: CSV file not found: {csv_file}")
        return
    
    print(f"Processing CSV file: {csv_file}")
    
    all_data = []
    
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)  
        print(f"CSV header: {header}")
        
        row_count = 0
        for row in reader:
            if len(row) < 4:
                print(f"Skipping incomplete row: {row}")
                continue
            
            base_name = row[0].strip()
            total_energy_str = row[1].strip()
            energy_per_atom_str = row[2].strip()
            cif_path = row[3].strip()
            
            try:
                total_energy = float(total_energy_str)
                energy_per_atom = float(energy_per_atom_str)
                
                chemical_formula = "Unknown"
                try:
                    structure = Structure.from_file(cif_path)
                    chemical_formula = structure.composition.formula.replace(" ", "")
                    print(f"Extracted chemical formula: {chemical_formula} from {cif_path}")
                except Exception as e:
                    print(f"Error reading CIF file {cif_path}: {e}")
                
                space_group = chem_dir_name
                
                all_data.append({
                    'chemical_formula': chemical_formula,
                    'space_group': space_group,
                    'base_name': base_name,
                    'total_energy': total_energy,
                    'energy_per_atom': energy_per_atom,
                    'cif_path': cif_path
                })
                row_count += 1
                
            except ValueError as e:
                print(f"Value conversion error: {e}, row data: {row}")
                continue
        
        print(f"Processed {row_count} rows from CSV file")
    
    if not all_data:
        print("No valid data found. Please check CSV file content and format.")
        return
    
    result_file_path = os.path.join(current_dir, '0_final_result_tot.txt')
    with open(result_file_path, 'w', newline='') as result_file:
        writer = csv.writer(result_file)
        writer.writerow(['Chemical_Formula', 'CIF_Base_Name', 'Total_Energy_eV', 'Energy_per_Atom_eV', 'Relaxed_CIF_Path', 'SG_ori'])
        
        for data in all_data:
            writer.writerow([
                data['chemical_formula'],
                data['base_name'],
                data['total_energy'],
                data['energy_per_atom'],
                data['cif_path'],
                data['space_group']
            ])
    
    print(f"\nProcessing completed:")
    print(f"Total structures processed: {len(all_data)}")
    print(f"Final result saved to: {result_file_path}")

if __name__ == "__main__":
    main()
