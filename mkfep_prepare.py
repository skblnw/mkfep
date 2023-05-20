import os
import sys
import subprocess
import shutil
import re
import numpy as np

def create_directories(rpath, positions, residues, structure_file, num_dirs, run_directory_name, use_predefined_replacements=True, check_directory_existence=True, free=True, complex=True):
    # Loop over all positions and residues
    for pos in positions:
        for res in residues:
            # Building the directory path string
            dir_path = os.path.join(rpath, f"pos{pos}/{res}")
            pdb2gmx_dir = os.path.join(dir_path, "pdb2gmx")
            run_dir = os.path.join(dir_path, run_directory_name)
            
            # Check if the directory already exists
            if check_directory_existence:
                if os.path.exists(dir_path):
                    # Ask the user if they want to continue if the directory already exists
                    response = input(f"The directory {dir_path} already exists. Do you want to continue? (y/n): ")
                    
                    # If the response is not 'y', skip this directory
                    if response.lower() != 'y':
                        print(f"Skipping directory: {dir_path}")
                        continue
                    else:
                        # Ask the user if they want to delete the existing directory
                        response = input(f"Do you want to delete the existing directory {dir_path}? (y/n): ")

                        # If the response is 'y', try to remove the existing directory and recreate it
                        # If any error occurs during this process, print the error message and skip this directory
                        if response.lower() == 'y':
                            try:
                                shutil.rmtree(dir_path)
                                os.makedirs(dir_path, exist_ok=True)
                            except Exception as e:
                                print(f"Error while removing and recreating directory {dir_path}: {e}")
                                continue
            
                # Check if pdb2gmx subdirectory already exists
                if not os.path.exists(pdb2gmx_dir):
                    # If it doesn't exist, create the pdb2gmx directories and run the pdb2fep script
                    create_pdb2gmx_directories_and_pdb2fep(rpath, pos, res, structure_file)
            
            # Check if run_directory_name already exists
            if os.path.exists(run_dir):
                # Ask the user if they want to continue if the directory already exists
                response = input(f"The directory {run_dir} already exists. Do you want to continue? (y/n): ")

                # If the response is not 'y', skip this directory
                if response.lower() != 'y':
                    print(f"Skipping directory: {run_dir}")
                    continue
            else:
                # Create the windows directories
                create_windows(rpath, pos, res, num_dirs, run_directory_name, use_predefined_replacements, free, complex)

def create_pdb2gmx_directories_and_pdb2fep(rpath, pos, res, structure_file):
    # Create the pdb2gmx directories
    # It is assumed that the 'rpath' directory contains directories named 'pos<pos>/<res>/pdb2gmx',
    # where <pos> and <res> are replaced with the function arguments
    pdb2gmx_dir = os.path.join(rpath, f"pos{pos}/{res}/pdb2gmx")
    os.makedirs(pdb2gmx_dir, exist_ok=True)

    # Write the position and residue information to the "mut" file
    with open(os.path.join(pdb2gmx_dir, "mut"), "w") as mut_file:
        mut_file.write(f"{pos} {res}")

    # Create symbolic links for the specified files in the pdb2gmx directory
    files_to_symlink = ['pdb2fep.py', structure_file]
    for file in files_to_symlink:
        src_path = os.path.join(rpath, file)
        os.symlink(src_path, os.path.join(pdb2gmx_dir, file))

    # Change the directory and run the pdb2fep script
    original_dir = rpath
    target_dir = pdb2gmx_dir  # reusing the previously defined path

    # The script that will be run in the target directory
    script = 'pdb2fep.py'

    # Define the log file path, the log file will be created in the target directory
    log_file_path = os.path.join(target_dir, 'pdb2fep.log')

    try:
        # Changing the working directory to the target directory
        os.chdir(target_dir)

        # Redirect stdout and stderr to the log file
        # This opens the log file in write mode, which means any existing file with the same name will be overwritten
        # If you want to append to an existing log file instead of overwriting, consider using mode 'a' instead of 'w'
        with open(log_file_path, 'w') as log_file:
            # Running the script using subprocess.run
            # The first argument to subprocess.run is a list where the first item is the command and the rest are arguments to the command
            # Here, "python" is the command and script is the argument
            # stdout is redirected to the log file
            # stderr is also redirected to the same log file by using subprocess.STDOUT
            # text=True means that the input and output are opened as text files
            subprocess.run(["python", script], stdout=log_file, stderr=subprocess.STDOUT, text=True)
    finally:
        # Changing the working directory back to the original directory
        # This is done in a finally block to ensure that the directory is changed back even if an error occurs
        os.chdir(original_dir)

def create_pdb2gmx_directories(rpath, pos, res, structure_file):
    # Build the directory path for pdb2gmx
    pdb2gmx_dir = os.path.join(rpath, f"pos{pos}/{res}/pdb2gmx")

    # Create the directory for pdb2gmx, 'exist_ok=True' allows the command to pass if the directory already exists
    os.makedirs(pdb2gmx_dir, exist_ok=True)

    # Open a file named "mut" in write mode in the pdb2gmx directory
    # The 'with' keyword is used here to handle the file. It automatically closes the file after the nested block of code
    with open(os.path.join(pdb2gmx_dir, "mut"), "w") as mut_file:
        # Write the position and residue to the "mut" file
        mut_file.write(f"{pos} {res}")

    # Define the list of files to create symbolic links for
    files_to_symlink = ['pdb2fep.py', structure_file]

    # Loop over each file in the files_to_symlink list
    for file in files_to_symlink:
        # Build the source path for the file
        src_path = os.path.join(rpath, file)

        # Create a symbolic link for the file in the pdb2gmx directory
        os.symlink(src_path, os.path.join(pdb2gmx_dir, file))

def change_directory_and_pdb2fep(rpath, pos, res):
    # Saving the original directory path
    original_dir = rpath

    # Building the target directory path string
    # It is assumed that the 'rpath' directory contains directories named 'pos<pos>/<res>/pdb2gmx',
    # where <pos> and <res> are replaced with the function arguments
    target_dir = os.path.join(rpath, f"pos{pos}/{res}/pdb2gmx")

    # The script that will be run in the target directory
    script = 'pdb2fep.py'

    # Define the log file path, the log file will be created in the target directory
    log_file_path = os.path.join(target_dir, 'pdb2fep.log')

    try:
        # Changing the working directory to the target directory
        os.chdir(target_dir)

        # Redirect stdout and stderr to the log file
        # This opens the log file in write mode, which means any existing file with the same name will be overwritten
        # If you want to append to an existing log file instead of overwriting, consider using mode 'a' instead of 'w'
        with open(log_file_path, 'w') as log_file:
            # Running the script using subprocess.run
            # The first argument to subprocess.run is a list where the first item is the command and the rest are arguments to the command
            # Here, "python" is the command and script is the argument
            # stdout is redirected to the log file
            # stderr is also redirected to the same log file by using subprocess.STDOUT
            # text=True means that the input and output are opened as text files
            subprocess.run(["python", script], stdout=log_file, stderr=subprocess.STDOUT, text=True)
    finally:
        # Changing the working directory back to the original directory
        # This is done in a finally block to ensure that the directory is changed back even if an error occurs
        os.chdir(original_dir)

def create_windows(rpath, position, residue, num_dirs, run_directory_name, use_predefined_replacements=True, free=True, complex=True):
    
    # Nested function to create a symlink for pdb2gmx in the target directory
    def symlink_pdb2gmx(original_dir, target_dir):

        pdb2gmx_source = os.path.relpath(os.path.join(original_dir, "pdb2gmx"), target_dir)
        pdb2gmx_target = os.path.join(target_dir, "pdb2gmx")

        # If the target file exists, remove it
        if os.path.exists(pdb2gmx_target):
            os.remove(pdb2gmx_target)

        # Create the symlink
        os.symlink(pdb2gmx_source, pdb2gmx_target)

    # Construct paths to original free and complex directories under pdb2gmx
    original_free_dir = os.path.join(rpath, f"pos{position}/{residue}/pdb2gmx/free")
    original_complex_dir = os.path.join(rpath, f"pos{position}/{residue}/pdb2gmx/complex")

    # Add free and complex options if respective flags are True
    options = []
    if free:
        options.append(('free', original_free_dir))
    if complex:
        options.append(('complex', original_complex_dir))

    # Loop over each option
    for option, original_dir in options:
        
        # Construct target directory path and create target directory if it doesn't exist
        target_dir = os.path.join(rpath, f"pos{position}/{residue}/{run_directory_name}/{option}")
        os.makedirs(target_dir, exist_ok=True)

        # Create symlink for pdb2gmx in target directory
        symlink_pdb2gmx(original_dir, target_dir)
        # Copy the 'mdp' directory into the target directory
        shutil.copytree(f"{rpath}/mdp", f"{rpath}/pos{position}/{residue}/{run_directory_name}/{option}", dirs_exist_ok=True)
        create_mdp_files(num_dirs, rpath, os.path.join(rpath, f"pos{position}/{residue}/{run_directory_name}/{option}"), use_predefined_replacements)
        create_run_files(rpath, position, residue, run_directory_name, option)

def create_mdp_files(num_dirs, rpath, base_dir, use_predefined_replacements=True):
    """
    Create MDP files in multiple directories with updated content.

    Args:
    num_dirs (int): Number of directories to create.
    rpath (str): Path to the original MDP files.
    base_dir (str): Path to the base directory where new directories will be created.
    """

    def get_replacement_strings(use_predefined_replacements=False, predefined_strings=None):
        if use_predefined_replacements and predefined_strings:
            init_lambda_state_replacement = predefined_strings.get('init_lambda_state')
            fep_lambdas_replacement = predefined_strings.get('fep_lambdas')
        else:
            total_number = predefined_strings.get('total_number')
            power = 1.1 # Adjust this value to control the intervals' distribution
            lambda_values = generate_symmetric_lambda_values(total_number, power)
            init_lambda_state_replacement = format_init_lambda_state_replacement(total_number)
            fep_lambdas_replacement = format_fep_lambdas_replacement(lambda_values)
        return init_lambda_state_replacement, fep_lambdas_replacement

    # Define the replacement strings
    # init_lambda_state_replacement = '; init_lambda_state        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23'
    # fep_lambdas_replacement = 'fep-lambdas = 0 0.00001 0.0001 0.001 0.01 0.02 0.04 0.06 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.94 0.96 0.98 0.99 0.999 0.9999 0.99999 1.00'
    # init_lambda_state_replacement = '; init_lambda_state        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47'
    # fep_lambdas_replacement = 'fep-lambdas              = 0 0.000001 0.00001 0.0001 0.001 0.004 0.01 0.02 0.03 0.04 0.05 0.06 0.08 0.1 0.12 0.14 0.18 0.22 0.26 0.30 0.34 0.38 0.42 0.46 0.5 0.54 0.58 0.62 0.66 0.7 0.74 0.78 0.82 0.86 0.88 0.9 0.92 0.94 0.95 0.96 0.97 0.98 0.99 0.996 0.999 0.9999 0.99999 0.999999 1.00'
    predefined_strings = {
        'total_number': num_dirs,
        'init_lambda_state': '; init_lambda_state        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27',
        'fep_lambdas': 'fep-lambdas = 0 0.00001 0.0001 0.001 0.01 0.02 0.03 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.92 0.94 0.96 0.97 0.98 0.99 0.999 0.9999 0.99999 1.00'
    }

    init_lambda_state_replacement, fep_lambdas_replacement = get_replacement_strings(use_predefined_replacements, predefined_strings=predefined_strings)

    for i in range(num_dirs):
        new_dir_path = os.path.join(base_dir, f"dir{i}")
        os.makedirs(new_dir_path, exist_ok=True)

        for prefix in ['em', 'nvt', 'md']:
            original_mdp_file_path = os.path.join(rpath, f"mdp/{prefix}.mdp")
            new_mdp_file_path = os.path.join(new_dir_path, f"{prefix}.mdp")

            with open(original_mdp_file_path, 'r') as original_mdp_file:
                mdp_content = original_mdp_file.read()

            mdp_content = re.sub(r'^; init_lambda_state.*', init_lambda_state_replacement, mdp_content, flags=re.MULTILINE)
            mdp_content = re.sub(r'^fep-lambdas.*', fep_lambdas_replacement, mdp_content, flags=re.MULTILINE)
            mdp_content = re.sub(r'^init_lambda_state.*', f'init_lambda_state = {i}', mdp_content, flags=re.MULTILINE)

            with open(new_mdp_file_path, 'w') as new_mdp_file:
                new_mdp_file.write(mdp_content)

def generate_symmetric_lambda_values(total_number, power):
    # Generate half lambda values first
    half_lambda_values = generate_half_lambda_values(total_number, power)
    mirrored_values = 1 - half_lambda_values[::-1]

    # Combine the two halves, excluding the middle value if the total number is odd
    if total_number % 2 == 0:
        lambda_values = np.concatenate((half_lambda_values, mirrored_values))
    else:
        lambda_values = np.concatenate((half_lambda_values, mirrored_values[1:]))
    return lambda_values

def generate_half_lambda_values(total_number, power):
    # Generate half of the total number of equally spaced values between 0 and 1
    half_total_number = (total_number + 1) // 2
    linear_values = np.linspace(0, .5, half_total_number)

    # Raise the linear values to the specified power
    scaled_values = linear_values ** power

    # Apply a non-linear transformation to obtain the desired distribution
    transformed_values = np.sin((scaled_values - 0.5) * np.pi) * 0.5 + 0.5
    return transformed_values

def format_init_lambda_state_replacement(total_number):
    indices = ' '.join(str(i) for i in range(total_number + 1))
    return f"; init_lambda_state {indices}"

def format_fep_lambdas_replacement(lambda_values):
    formatted_values = ' '.join(f"{value:.5f}" for value in lambda_values)
    return f"fep-lambdas = {formatted_values}"

def print_values_and_intervals(lambda_values):
    print("Lambda values:")
    for value in lambda_values:
        print(f"{value:.4f}")
    intervals = np.diff(lambda_values)
    for i, interval in enumerate(intervals):
        print(f"Interval {i + 1}: {interval:.4f}")

def create_pbs_files(rpath, pos, res, run_directory_name, option):
    with open(os.path.join(rpath, "pbs"), "r") as pbs_file:
        pbs_content = pbs_file.read()
    pbs_content = re.sub(r'#PBS -N.*', f"#PBS -N p{pos}-{res}-{option}", pbs_content)
    with open(os.path.join(rpath, f"pos{pos}/{res}/{run_directory_name}/{option}/pbs"), "w") as pbs_file:
        pbs_file.write(pbs_content)
    # subprocess.run(["qsub", os.path.join(rpath, f"pos{pos}/{res}/{run_directory_name}/{option}/pbs")])

def create_run_files(rpath, pos, res, run_directory_name, option):
    with open(os.path.join(rpath, "run"), "r") as run_file:
        run_content = run_file.read()
    with open(os.path.join(rpath, f"pos{pos}/{res}/{run_directory_name}/{option}/run"), "w") as run_file:
        run_file.write(run_content)
    # subprocess.run(["bash", os.path.join(rpath, f"pos{pos}/{res}/{run_directory_name}/{option}/run")])

def structure_files_exist(structure_files=None):
    if structure_files:
        for structure_file in structure_files:
            if os.path.isfile(structure_file):
                return structure_file

    print("Error: Neither 'md.gro' nor 'md.pdb' were found.")
    sys.exit(1)

def check_existence(required_files=None, required_directories=None):
    if required_files:
        for file in required_files:
            if not os.path.isfile(file):
                print(f"Error: Required file '{file}' not found.")
                sys.exit(1)

    if required_directories:
        for directory in required_directories:
            if not os.path.isdir(directory):
                print(f"Error: Required directory '{directory}' not found.")
                sys.exit(1)

def check_command_existence(command_name):
    if not shutil.which(command_name):
        print(f"Error: {command_name} does not exist", file=sys.stderr)
        sys.exit(1)

def one_to_three(aa_seq):
    aa_dict = {
        "A": "ALA",  # Alanine
        "R": "ARG",  # Arginine
        "N": "ASN",  # Asparagine
        "D": "ASP",  # Aspartic acid
        "C": "CYS",  # Cysteine
        "Q": "GLN",  # Glutamine
        "E": "GLU",  # Glutamic acid
        "G": "GLY",  # Glycine
        "H": "HIS",  # Histidine
        "I": "ILE",  # Isoleucine
        "L": "LEU",  # Leucine
        "K": "LYS",  # Lysine
        "M": "MET",  # Methionine
        "F": "PHE",  # Phenylalanine
        "P": "PRO",  # Proline
        "S": "SER",  # Serine
        "T": "THR",  # Threonine
        "W": "TRP",  # Tryptophan
        "Y": "TYR",  # Tyrosine
        "V": "VAL"   # Valine
    }
    # Use dictionary to map each one-letter code to its corresponding three-letter code
    # Apply mapping to each character in the input string using a list comprehension
    # Call the upper() method on the input string to ensure that all characters are in uppercase
    return [aa_dict[aa] for aa in aa_seq.upper()]

def main():
    rpath = os.getcwd()
    positions = [2]
    residues = one_to_three("A")
    # residues = one_to_three("ARNDCQEHILKMFSTWYV")
    num_windows = 10
    run_directory_name = "win10.t2"
    pdb2gmx_needed = False
    use_predefined_replacements = False

    # Check if required files exist
    structure_file = structure_files_exist(["md.gro", "md.pdb"])

    if pdb2gmx_needed:
        # Check if required command exist
        check_command_existence("pmx")

        # Check if required files exist
        check_existence(
            required_files=["pdb2fep.py", "run"],
            required_directories=["mdp"]
        )

        create_directories(rpath, positions, residues, structure_file, num_windows, run_directory_name, use_predefined_replacements, check_directory_existence=pdb2gmx_needed, free=True, complex=True)
    else:
        print("You specified `no pdb2gmx needed`. Skipping directory creation and file existence check.")

        create_directories(rpath, positions, residues, structure_file, num_windows, run_directory_name, use_predefined_replacements, check_directory_existence=pdb2gmx_needed, free=True, complex=True)

if __name__ == "__main__":

    main()
