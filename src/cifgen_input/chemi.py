import math
from functools import reduce

def gcd_multiple(*numbers):
    return reduce(math.gcd, numbers)

def simplify_formula(*values):
    gcd_value = gcd_multiple(*values)
    return [value // gcd_value for value in values]

def parse_value_spec(value_str):
    """
    Parse value specification string.
    Formats:
      - "6"   -> range mode: 1-6 (backward compatible)
      - "=6"  -> fixed mode: only 6
      - "2-6" -> range mode: 2-6
    Returns: list of values
    """
    value_str = value_str.strip()
    
    if value_str.startswith('='):
        fixed_val = int(value_str[1:])
        return [fixed_val]
    elif '-' in value_str:
        parts = value_str.split('-')
        start_val = int(parts[0])
        end_val = int(parts[1])
        return list(range(start_val, end_val + 1))
    else:
        max_value = int(value_str)
        return list(range(1, max_value + 1))

def generate_chemical_formulas(elements_list):
    formulas = set()
    
    def helper(index, current_elements, current_values):
        if index == len(elements_list):
            simplified_values = simplify_formula(*current_values)
            formula = ''.join(f"{element}{value}" for element, value in zip(current_elements, simplified_values))
            formulas.add(formula)
            return
        
        elements, values = elements_list[index]
        for element in elements:
            for value in values:
                current_elements.append(element)
                current_values.append(value)
                helper(index + 1, current_elements, current_values)
                current_elements.pop()
                current_values.pop()
    
    helper(0, [], [])
    return sorted(formulas)

def save_to_file(formulas, filename):
    with open(filename, 'w') as file:
        for formula in formulas:
            file.write(formula + '\n')

import os

def read_chem_input(filename):
    if os.path.exists(filename):
        with open(filename, 'r') as file:
            lines = file.readlines()
    else:
        base_name = filename
        extensions = ['.txt', '.inp', '.in', '', '.dat']
        
        found = False
        for ext in extensions:
            try_name = base_name + ext if not base_name.endswith(ext) else base_name
            if os.path.exists(try_name):
                with open(try_name, 'r') as file:
                    lines = file.readlines()
                found = True
                break
        
        if not found:
            raise FileNotFoundError(f"Cannot find chem_input file with any supported extension: {extensions}")
    
    elements_list = []
    for i in range(0, len(lines), 2):
        elements = lines[i].strip().split(',')
        value_spec = lines[i + 1].strip()
        values = parse_value_spec(value_spec)
        elements_list.append((elements, values))
    
    return elements_list
  
def main(): 
    elements_list = read_chem_input('chem_input')
    formulas = generate_chemical_formulas(elements_list)
    save_to_file(formulas, 'oxidation')

if __name__ == "__main__":  
    main()