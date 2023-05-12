import pandas as pd
import numpy as np
import sys

# Read data from the file and create a DataFrame
def read_data(file_path):
    data_list = []

    with open(file_path, "r") as file:
        for line in file:
            if "broken" in line:
                continue
            if "not exist" in line:
                continue
            tokens = line.strip().split()
            position = int(tokens[0])
            aa = tokens[1].strip(":")
            ddg = float(tokens[2])
            se = float(tokens[4].strip("+-"))
            data_list.append({"Position": position, "AA": aa, "ddG": ddg, "se": se})

    return pd.DataFrame(data_list)

# Create a mapping of amino acid abbreviations to single-letter codes
def create_mapping():
    letter_list = list('ACDEFGHIKLMNPQRSTVWY')
    res_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
    return {res_list[ii]: letter for ii, letter in enumerate(letter_list)}

# Mutate the reference sequence at the given position with the target amino acid
def mutate_ref(ref, pos, target):
    ref_list = list(ref)
    ref_list[pos - 1] = target
    return ''.join(ref_list)

def main(file_path):
    data = read_data(file_path)

    ref0 = 'ITDQVPFSV'
    aa_mapping = create_mapping()

    data = data.dropna(subset=['ddG'])
    data['target'] = data['AA'].map(aa_mapping)
    data['target'] = data.apply(lambda row: mutate_ref(ref0, row['Position'], row['target']), axis=1)

    # Insert 'ref0' as the first column in the output
    data.insert(0, 'original', ref0)
    data[['original', 'target', 'ddG', 'se']].to_csv('ddg_formatted.csv', header=['original', 'target', 'ddg', 'ddgse'], index=False)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("mkpp> Please provide the path to the data file (e.g., 'combined_results.txt')")
    else:
        file_path = sys.argv[1]
        main(file_path)
