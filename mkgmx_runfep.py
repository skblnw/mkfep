import os
import sys
import subprocess
import shutil
import re

def create_directories(rpath, positions, residues, structure_file, num_dirs, run_directory_name, free=True, complex=True):
    for pos in positions:
        for res in residues:
            dir_path = os.path.join(rpath, f"pos{pos}/{res}")
            
            # Check if the directory already exists
            if os.path.exists(dir_path):
                response = input(f"The directory {dir_path} already exists. Do you want to continue? (y/n): ")
                
                if response.lower() != 'y':
                    print(f"Skipping directory: {dir_path}")
                    continue
                else:
                    try:
                        shutil.rmtree(dir_path)
                        os.makedirs(dir_path, exist_ok=True)
                    except Exception as e:
                        print(f"Error while removing and recreating directory {dir_path}: {e}")
                        continue
            create_pdb2gmx_directories(rpath, pos, res, structure_file)
            create_windows(rpath, pos, res, num_dirs, run_directory_name, free, complex)
            change_directory_and_run_script(rpath, pos, res)

def create_pdb2gmx_directories(rpath, pos, res, structure_file):
    pdb2gmx_dir = os.path.join(rpath, f"pos{pos}/{res}/pdb2gmx")
    os.makedirs(pdb2gmx_dir, exist_ok=True)

    with open(os.path.join(pdb2gmx_dir, "mut"), "w") as mut_file:
        mut_file.write(f"{pos} {res}")

    src_path = os.path.join(os.path.expanduser("~"), "charmm36-feb2021.ff")
    os.symlink(src_path, os.path.join(pdb2gmx_dir, "charmm36-feb2021.ff"))

    files_to_symlink = ['mkgmx_pdb2fep.py', structure_file]
    for file in files_to_symlink:
        src_path = os.path.join(rpath, file)
        os.symlink(src_path, os.path.join(pdb2gmx_dir, file))

def create_windows(rpath, position, residue, num_dirs, run_directory_name, free=True, complex=True):
    def symlink_pdb2gmx(original_dir, target_dir):
        """Create a symlink for pdb2gmx in the target directory."""
        pdb2gmx_source = os.path.relpath(os.path.join(original_dir, "pdb2gmx"), target_dir)
        pdb2gmx_target = os.path.join(target_dir, "pdb2gmx")
        os.symlink(pdb2gmx_source, pdb2gmx_target)

    original_free_dir = os.path.join(rpath, f"pos{position}/{residue}/pdb2gmx/free")
    original_complex_dir = os.path.join(rpath, f"pos{position}/{residue}/pdb2gmx/complex")

    options = []
    if free:
        options.append(('free', original_free_dir))
    if complex:
        options.append(('complex', original_complex_dir))

    for option, original_dir in options:
        target_dir = os.path.join(rpath, f"pos{position}/{residue}/{run_directory_name}/{option}")
        os.makedirs(target_dir, exist_ok=True)

        symlink_pdb2gmx(original_dir, target_dir)
        shutil.copytree(f"{rpath}/mdp", f"{rpath}/pos{position}/{residue}/{run_directory_name}/{option}", dirs_exist_ok=True)
        create_mdp_files(num_dirs, rpath, os.path.join(rpath, f"pos{position}/{residue}/{run_directory_name}/{option}"))
        create_run_files(rpath, position, residue, run_directory_name, option)

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

def generate_half_lambda_values(total_number, power):
    # Generate half of the total number of equally spaced values between 0 and 1
    half_total_number = (total_number + 1) // 2
    linear_values = np.linspace(0, .5, half_total_number)

    # Raise the linear values to the specified power
    scaled_values = linear_values ** power

    # Apply a non-linear transformation to obtain the desired distribution
    transformed_values = np.sin((scaled_values - 0.5) * np.pi) * 0.5 + 0.5

    return transformed_values

def generate_symmetric_lambda_values(total_number, power):
    half_lambda_values = generate_half_lambda_values(total_number, power)
    mirrored_values = 1 - half_lambda_values[::-1]

    # Combine the two halves, excluding the middle value if the total number is odd
    if total_number % 2 == 0:
        lambda_values = np.concatenate((half_lambda_values, mirrored_values))
    else:
        lambda_values = np.concatenate((half_lambda_values, mirrored_values[1:]))

    return lambda_values

def print_values_and_intervals(lambda_values):
    print("Lambda values:")
    for value in lambda_values:
        print(f"{value:.4f}")
    intervals = np.diff(lambda_values)
    for i, interval in enumerate(intervals):
        print(f"Interval {i + 1}: {interval:.4f}")

def format_init_lambda_state_replacement(total_number):
    indices = ' '.join(str(i) for i in range(total_number + 1))
    return f"; init_lambda_state {indices}"

def format_fep_lambdas_replacement(lambda_values):
    formatted_values = ' '.join(f"{value:.5f}" for value in lambda_values)
    return f"fep-lambdas = {formatted_values}"

def get_replacement_strings(use_predefined_strings=False, predefined_strings=None):
    if use_predefined_strings and predefined_strings:
        init_lambda_state_replacement = predefined_strings.get('init_lambda_state')
        fep_lambdas_replacement = predefined_strings.get('fep_lambdas')
    else:
        total_number = predefined_strings.get('total_number')
        power = 1.1 # Adjust this value to control the intervals' distribution
        lambda_values = generate_symmetric_lambda_values(total_number, power)
        init_lambda_state_replacement = format_init_lambda_state_replacement(total_number)
        fep_lambdas_replacement = format_fep_lambdas_replacement(lambda_values)
    
    return init_lambda_state_replacement, fep_lambdas_replacement

def create_mdp_files(num_dirs, rpath, base_dir):
    """
    Create MDP files in multiple directories with updated content.

    Args:
    num_dirs (int): Number of directories to create.
    rpath (str): Path to the original MDP files.
    base_dir (str): Path to the base directory where new directories will be created.
    """

    # Define the replacement strings
    # init_lambda_state_replacement = '; init_lambda_state        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23'
    # fep_lambdas_replacement = 'fep-lambdas = 0 0.00001 0.0001 0.001 0.01 0.02 0.04 0.06 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.94 0.96 0.98 0.99 0.999 0.9999 0.99999 1.00'
    predefined_strings = {
        'total_number': num_dirs,
        'init_lambda_state': '; init_lambda_state        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27',
        'fep_lambdas': 'fep-lambdas = 0 0.00001 0.0001 0.001 0.01 0.02 0.03 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.92 0.94 0.96 0.97 0.98 0.99 0.999 0.9999 0.99999 1.00'
    }

    init_lambda_state_replacement, fep_lambdas_replacement = get_replacement_strings(use_predefined_strings=True, predefined_strings=predefined_strings)

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

def change_directory_and_run_script(rpath, pos, res):
    original_dir = rpath
    target_dir = os.path.join(rpath, f"pos{pos}/{res}/pdb2gmx")

    script = 'mkgmx_pdb2fep.py'
    # Define the log file path
    log_file_path = os.path.join(target_dir, 'pdb2fep.log')

    try:
        os.chdir(target_dir)
        # Redirect stdout and stderr to the log file
        with open(log_file_path, 'w') as log_file:
            subprocess.run(["python", script], stdout=log_file, stderr=subprocess.STDOUT, text=True)
    finally:
        os.chdir(original_dir)

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

def main(structure_file):
    rpath = os.getcwd()
    positions = [2,3,4,5,6,7,8]
    residues = one_to_three("A")
    num_windows = 28
    run_directory_name = "win28.t1"
    create_directories(rpath, positions, residues, structure_file, num_windows, run_directory_name, free=True, complex=True)

if __name__ == "__main__":

    # Check if required command exist
    check_command_existence("pmx")

    # Check if required files exist
    structure_file = structure_files_exist(["md.gro", "md.pdb"])
    check_existence(
        required_files=["mkgmx_pdb2fep.py", "run"],
        required_directories=["mdp"]
    )

    main(structure_file)
