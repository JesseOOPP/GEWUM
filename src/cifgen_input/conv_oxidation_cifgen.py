import os

ELEMENTS = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
    'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs',
    'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
]

def read_possible_oxidation_states(file_path):
    compounds_with_oxidation = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        i = 0
        while i < len(lines):
            if lines[i].startswith('Compound:'):
                compound = lines[i].strip().split(': ')[1].replace(' ', '')  
                i += 1
                if i < len(lines) and lines[i].strip().startswith('Possible Oxidation States:'):
                    i += 1
                    if i < len(lines) and lines[i].strip().startswith('{'):
                        compounds_with_oxidation.append(compound)
                        print(f"Added compound: {compound}") 
            else:
                i += 1
    return compounds_with_oxidation

def filter_oxidation_inp(input_file, output_file, compounds_with_oxidation):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            compound = line.strip()
            if compound in compounds_with_oxidation:
                outfile.write(line)
                print(f"Kept compound: {compound}")

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

def cleanup_files(files_to_delete):
    for file_path in files_to_delete:
        if os.path.exists(file_path):
            os.remove(file_path)
            print(f"Deleted {file_path}")

def main():
    possible_oxidation_file = 'possible_oxidation_states.txt'
    input_file = 'oxidation'
    filtered_oxidation_file = 'filtered_oxidation'
    cifgen_file = 'cifgen.inp'
    
    files_to_delete = [input_file, filtered_oxidation_file]
    
    try:
        compounds_with_oxidation = read_possible_oxidation_states(possible_oxidation_file)
        print(f"Compounds with oxidation states: {compounds_with_oxidation}")
        
        filter_oxidation_inp(input_file, filtered_oxidation_file, compounds_with_oxidation)
        print(f"Filtered compounds have been written to {filtered_oxidation_file}")
        
        compounds = read_oxidation_file(filtered_oxidation_file)
        
        user_input = input("Input a max number in cif generation per SG:")
        
        cifgen_lines = convert_to_cifgen_format(compounds, user_input)
        
        write_cifgen_file(cifgen_file, cifgen_lines)
        
    finally:
        cleanup_files(files_to_delete)

if __name__ == "__main__":
    main()
