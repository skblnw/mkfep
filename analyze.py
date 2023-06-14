"""
# Topology File Analyzer

This script analyzes topology files typically used in molecular dynamics simulations. The files are expected to have multiple sections such as `moleculetype`, `atoms`, `bonds`, etc., each starting with a section header enclosed in square brackets (`[ ]`).

The script performs the following tasks:

1. Reads a given topology file line by line.
2. Identifies and prints the total number of sections and their names.
3. Within each section, counts the number of entries (lines) and groups them based on the number of columns. Comment lines starting with `;` or `#` are ignored.
4. It treats the group with the smallest number of columns as the 'original' group, and all other groups as 'perturbed'.
5. Writes the lines belonging to the 'original' and 'perturbed' groups into separate output files, preserving the section headers.

## Usage

```shell
python script_name.py --file_path <input_file> --original_file_name <original_output_file> --perturbed_file_name <perturbed_output_file>
```
The --file_path argument is mandatory and should be the path to the input topology file. 
The --original_file_name and --perturbed_file_name arguments are optional; if not provided, the script will write to 'original.txt' and 'perturbed.txt' by default.

##Functionality
The main function analyze_file(file_path, original_file_name, perturbed_file_name) performs the tasks mentioned above. 
It opens the input file and reads it line by line. 
When it encounters a section header, it stores the section name and begins counting the number of entries and columns for that section. 
The lines are temporarily stored in memory until the script encounters a new section or reaches the end of the file. 
At this point, the function identifies the group with the smallest number of columns ('original') and writes the lines to the corresponding output files. 
The process is repeated for all sections in the file.
"""

import sys
import argparse
from collections import defaultdict

def analyze_file(file_path, original_file_name, perturbed_file_name):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    section_names = []
    current_section = None
    column_counts = defaultdict(int)
    section_lines = defaultdict(list)
    min_columns = float('inf')
    max_columns = float('-inf')
    original_file = open(original_file_name, 'w')
    perturbed_file = open(perturbed_file_name, 'w')

    for line in lines:
        stripped_line = line.strip()
        if stripped_line.startswith('[') and stripped_line.endswith(']'):
            if current_section is not None:
                # Print column counts for each section
                print(f"\n{current_section} column counts: ")
                for columns, count in column_counts.items():
                    print(f"{columns} columns: {count} entries")
                # Write lines to appropriate files
                for section_line in section_lines[min_columns]:
                    original_file.write(section_line)
                for columns in range(min_columns+1, max_columns+1):
                    for section_line in section_lines[columns]:
                        perturbed_file.write(section_line)
                section_names.append(current_section)
                # Reset for new section
                column_counts = defaultdict(int)
                section_lines = defaultdict(list)
                min_columns = float('inf')
                max_columns = float('-inf')

            current_section = stripped_line[1:-1].strip()  # Remove brackets and strip whitespace
            original_file.write(line)
            perturbed_file.write(line)
        elif not stripped_line.startswith(';') and not stripped_line.startswith('#') and stripped_line != '':
            columns = len(stripped_line.split())
            column_counts[columns] += 1
            section_lines[columns].append(line)
            min_columns = min(min_columns, columns)
            max_columns = max(max_columns, columns)

    # Process the last section
    if current_section is not None:
        # Print column counts for each section
        print(f"\n{current_section} column counts: ")
        for columns, count in column_counts.items():
            print(f"{columns} columns: {count} entries")
        # Write lines to appropriate files
        for section_line in section_lines[min_columns]:
            original_file.write(section_line)
        for columns in range(min_columns+1, max_columns+1):
            for section_line in section_lines[columns]:
                perturbed_file.write(section_line)
        section_names.append(current_section)

    print(f'\nThere are {len(section_names)} sections in the file.')
    
    original_file.close()
    perturbed_file.close()

# Argument parsing
parser = argparse.ArgumentParser(description='Analyze a topology file.')
parser.add_argument('file_path', type=str, help='Path to the topology file to analyze.')
parser.add_argument('--original_file_name', type=str, default='original.txt', help='Name of the output file for original lines.')
parser.add_argument('--perturbed_file_name', type=str, default='perturbed.txt', help='Name of the output file for perturbed lines.')
args = parser.parse_args()

#Use the function on your file
analyze_file(args.file_path, args.original_file_name, args.perturbed_file_name)
