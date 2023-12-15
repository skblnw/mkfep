"""
# Python Script for Extracting Amino Acid Sequence from a PDB File

This Python script extracts the amino acid sequence of a specified chain from a PDB file.

**How It Works:**

1. **Define a Code Dictionary:** The script contains a dictionary mapping three-letter amino acid codes to one-letter codes.
2. **Parse the PDB File:** The `get_sequence_from_pdb` function reads the PDB file line by line. It selects lines starting with "ATOM" or "HETATM" and corresponding to the specified chain.
3. **Extract Amino Acid Codes:** For each selected line, it extracts the three-letter amino acid code and converts it to a one-letter code using the dictionary.
4. **Sequence Creation:** It forms the sequence by appending each one-letter code, ensuring each residue is added only once.
5. **Command Line Arguments:** The `main` function uses the argparse library to handle command line inputs: the PDB filename and the chain id.
6. **Print the Sequence:** After extracting the sequence, the script prints it out.

To use the script from the command line, the syntax is:
```
python get_sequence.py my_pdb_file.pdb A
```

This will output the sequence for chain A from `my_pdb_file.pdb`.

"""

import argparse

# Residue names in three-letter format and their corresponding one-letter code
AA_DICT = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def get_sequence_from_pdb(filename, chain):
    sequence = []
    with open(filename, 'r') as file:
        lines = file.readlines()
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            residue = line[17:20].strip()
            chain_id = line[21].strip()
            res_seq = int(line[22:26].strip())
            if chain_id == chain and residue in AA_DICT:
                if not sequence or sequence[-1][1] != res_seq:
                    sequence.append((AA_DICT[residue], res_seq))
    return "".join([residue[0] for residue in sequence])

def main():
    parser = argparse.ArgumentParser(description='Get the amino acid sequence of a specific chain from a PDB file.')
    parser.add_argument('filename', type=str, help='The pdb filename.')
    parser.add_argument('chain', type=str, help='The chain id.')
    args = parser.parse_args()

    sequence = get_sequence_from_pdb(args.filename, args.chain)
    print(f"The sequence for chain {args.chain} is: {sequence}")

if __name__ == "__main__":
    main()
