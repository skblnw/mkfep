#!/bin/python
# python edit_lambda_values.py xxx.fepout

import re, sys


def edit_lambda_values(filename, output_filename):
    """
    Edits the lambda values in an FEP output file to 1 - original value and rounds to 8 decimal places.
    This version also considers numbers in scientific notation.

    Args:
    filename (str): The path to the FEP output file.
    output_filename (str): The path to save the edited file.
    """
    number_pattern = re.compile(r"\d+(\.\d+)?(e-?\d+)?")  # Pattern to match numbers including scientific notation

    with open(filename, 'r') as file, open(output_filename, 'w') as outfile:
        for line in file:
            if "LAMBDA SET TO" in line or "Free energy change for lambda window" in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if number_pattern.match(part) and "STEPS" not in line:
                        lambda_val = float(part)
                        parts[i] = str(round(1-lambda_val, 8))
                outfile.write(' '.join(parts)+'\n')
            elif "EQUILIBRATION AT LAMBDA" in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if number_pattern.match(part):
                        lambda_val = float(part)
                        parts[i] = str(round(1-lambda_val, 8))
                outfile.write(' '.join(parts)+'\n')
            else:
                outfile.write(line)


if __name__ == '__main__':
    file = sys.argv[1]
    edit_lambda_values(file, file+'.fepout')

# test
# edit_lambda_values('/data/work/make_hybrid_top/tests/Lei2Guanqiao/D3N/bound/t1/alchemy.fepout', '/data/work/make_hybrid_top/tests/Lei2Guanqiao/D3N/bound/t1/alchemy.fepout.fepout')
