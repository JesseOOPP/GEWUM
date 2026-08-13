import math
import re


def calculate_gcd_for_formula(formula):
    pattern = r'([A-Z][a-z]*)(\d*)'
    counts = []
    for element, count_str in re.findall(pattern, formula):
        count = int(count_str) if count_str else 1
        counts.append(count)
    
    if not counts:
        return 1
    gcd_val = counts[0]
    for num in counts[1:]:
        gcd_val = math.gcd(gcd_val, num)
    return gcd_val


def main():
    data = []
    with open('energy_final.txt', 'r') as f:
        next(f)
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            if len(parts) < 6: 
                continue
            formula = parts[0].strip()  
            filename = parts[1].strip() 
            try:
                energy = float(parts[2].strip()) 
            except ValueError:
                continue
            relaxed_cif_path = parts[4].strip() 
            data.append((formula, filename, energy, relaxed_cif_path))

    processed_data = []
    for formula, filename, energy, relaxed_cif_path in data:
        gcd_val = calculate_gcd_for_formula(formula)
        adjusted_energy = energy / gcd_val
        processed_data.append((formula, filename, adjusted_energy, relaxed_cif_path))

    processed_data.sort(key=lambda x: x[2]) 

    with open('sorted_energy_final.txt', 'w') as f_out:
        f_out.write("Chemical_Formula,Relaxed_CIF_Path,Energy_eV_per_Atom\n")
        for formula, filename, adjusted_energy, relaxed_cif_path in processed_data:
            formatted_energy = f"{adjusted_energy:.3g}"
            f_out.write(f"{formula},{relaxed_cif_path},{formatted_energy}\n")
    
    print(f"Reorder complete: {len(processed_data)} entries written to sorted_energy_final.txt")


if __name__ == "__main__":
    main()
