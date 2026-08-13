import os
import re

properties = {
    "Bulk Modulus B (GPa)": "B",
    "Young's Modulus E (GPa)": "E",
    "Shear Modulus G (GPa)": "G",
    "Poisson's Ratio v": "V",
    "P-wave Modulus (GPa)": "PM",
    "Pugh's Ratio (B/G)": "B/G",
    "Vickers Hardness (GPa)[6]": "VH6",
    "Vickers Hardness (GPa)[7]": "VH7"
}


def main():
    output_file = "ela_tot.dat"
    with open(output_file, 'w') as out_f:
        headers = ["Directory"] + list(properties.values())
        out_f.write("\t".join(headers) + "\n")
        
        for dir_name in os.listdir('.'):
            if os.path.isdir(dir_name):
                ela_file = os.path.join(dir_name, "ela.dat")
                
                if os.path.exists(ela_file):
                    try:
                        with open(ela_file, 'r') as f:
                            content = f.read()
                        
                        avg_section_start = content.find("Average mechanical properties of bulk polycrystal:")                    
                        avg_section = content[avg_section_start:]
                        
                        hill_values = []
                        for prop in properties:
                            pattern = re.compile(rf"\| *{re.escape(prop)} *\|.*?\|.*?\|.*?\|.*?\|.*?\| *([\d.]+) *\|")
                            match = pattern.search(avg_section)
                            
                            if match:
                                hill_values.append(match.group(1))
                            else:
                                alt_pattern = re.compile(rf"\| *{re.escape(prop)} .*?\|.*?\|.*?\| *([\d.]+) *\|")
                                alt_match = alt_pattern.search(avg_section)
                                if alt_match:
                                    hill_values.append(alt_match.group(1))
                                else:
                                    hill_values.append("N/A")
                        
                        out_f.write(f"{dir_name}\t" + "\t".join(hill_values) + "\n")
                    
                    except Exception as e:
                        print(f"Error with {dir_name} : {str(e)}")
                else:
                    print(f"{dir_name} Warning: no ela.dat file")

    print(f"All info extracted {output_file}")


if __name__ == "__main__":
    main()
