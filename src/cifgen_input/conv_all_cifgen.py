import os

def read_oxidation_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return [line.strip() for line in lines]

def parse_compound(compound):
    elements = []
    counts = []
    i = 0
    while i < len(compound):
        if compound[i].isupper():
            if i + 1 < len(compound) and compound[i + 1].islower():
                element = compound[i:i+2]
                i += 2
            else:
                element = compound[i]
                i += 1
            elements.append(element)
            count = 0
            while i < len(compound) and compound[i].isdigit():
                count = count * 10 + int(compound[i])
                i += 1
            counts.append(count if count > 0 else 1)
        else:
            raise ValueError(f"Invalid character '{compound[i]}' at position {i} in compound '{compound}'")
    return elements, counts

def convert_to_cifgen_format(compounds, user_input):
    cifgen_lines = []
    for compound in compounds:
        elements, counts = parse_compound(compound)
        cifgen_line = f"{elements},{tuple(counts)},{user_input}"
        cifgen_lines.append(cifgen_line)
    return cifgen_lines

def write_cifgen_file(file_path, cifgen_lines):
    with open(file_path, 'w') as file:
        for line in cifgen_lines:
            file.write(line + '\n')

def main():
    input_file = 'oxidation'
    cifgen_file = 'cifgen.inp'
    
    try:
        compounds = read_oxidation_file(input_file)
        user_input = input("Input a max number in cif generation per Space Group:")
        cifgen_lines = convert_to_cifgen_format(compounds, user_input)
        write_cifgen_file(cifgen_file, cifgen_lines)
        
    finally:
        if os.path.exists(input_file):
            os.remove(input_file)

if __name__ == "__main__":
    main()