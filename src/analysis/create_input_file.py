# Searches all .fepout files. Comment those you don't want.
# Usage: python create_input_file.py -f forward_folder -b backward_folder
import os
import argparse
import glob


def create_input_file(forward, backward):
    fepout_file = f'{forward}-{backward}.csv'
    if os.path.exists(fepout_file):
        os.remove(fepout_file)
    with open(fepout_file, 'w') as file:
        file.write("# Comment lines with #\n")

    search_pattern = os.path.join(os.getcwd(), "**", "*.fepout")  # any fepout files here
    for filepath in sorted(glob.glob(search_pattern, recursive=True)):
        if 'equil' in filepath:
            continue
        else:
            normalized_path = os.path.normpath(filepath)  # Normalize path for Windows

        if "complex" in normalized_path or 'bound' in normalized_path:
            state = "bound"
        elif "ligand" in normalized_path or 'free' in normalized_path:
            state = "free"
        else:
            continue  # Skip files not in the desired directories

        if f"{forward}" in normalized_path:
            direction = "forward"
        elif f"{backward}" in normalized_path:
            direction = "backward"
        else:
            continue  # Skip files not matching forward or backward

        # Convert paths to Windows format for the output
        windows_path = filepath.replace(os.sep, '\\')
        with open(fepout_file, 'a') as file:
            file.write(f"{state},{direction},{windows_path}\n")

    print("Done, and run with:")
    print(f"python fepout_analysis.py -f {fepout_file}")
    print('Without proper prefixes or starting with #, lines will be skipped.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create input file for FEP analysis.')
    parser.add_argument('-f','--forward', type=str, help='Forward transformation identifier.')
    parser.add_argument('-b','--backward', type=str, help='Backward transformation identifier.')
    args = parser.parse_args()

    create_input_file(args.forward, args.backward)
