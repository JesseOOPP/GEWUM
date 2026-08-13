from pymatgen.core import Composition

def process_oxidation(input_file="oxidation", output_file="possible_oxidation_states.txt"):
    compounds = []
    with open(input_file, "r") as file:
        for line in file:
            compound_str = line.strip()
            if not compound_str:
                continue
            try:
                compounds.append(Composition(compound_str))
            except Exception as e:
                print(f"Error processing line '{line.strip()}': {e}")

    with open(output_file, "w") as file:
        for compound in compounds:
            oxi_states = compound.oxi_state_guesses(all_oxi_states=False)
            file.write(f"Compound: {compound}\nPossible Oxidation States:\n")
            for oxi_state in oxi_states:
                file.write(f"  {oxi_state}\n")
            file.write("\n")
    print(f"Results written to {output_file}")

def main(): 
    process_oxidation()

if __name__ == "__main__": 
    main()
