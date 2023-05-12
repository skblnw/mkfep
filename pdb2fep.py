import os
import sys
import subprocess
import shutil
import glob
from pathlib import Path

# Variables
PREFIX = "md"
OUTPUT_FILE = f"{PREFIX}_with_chains.pdb"
OUTPUT_DIR = "chains"

# Set INPUT_FILE
gro_file = f"{PREFIX}.gro"
pdb_file = f"{PREFIX}.pdb"
if os.path.exists(gro_file):
    INPUT_FILE = gro_file
elif os.path.exists(pdb_file):
    INPUT_FILE = pdb_file
else:
    raise FileNotFoundError(f"Neither {pdb_file} nor {gro_file} exist.")

# Check if required files exist
required_files = [
    "mut"
]
for file in required_files:
    if not os.path.isfile(file):
        print(f"Error: Required file '{file}' not found.")
        sys.exit(1)

# Convert gro to pdb if needed
if INPUT_FILE.endswith(".gro"):
    subprocess.run(["gmx", "editconf", "-f", INPUT_FILE, "-o", f"{PREFIX}.pdb"], check=True)
    INPUT_FILE = f"{PREFIX}.pdb"

def increment_chain(chain):
    return chr(ord(chain) + 1)

# Assign chain identifiers to repeating residue number ranges
with open(INPUT_FILE, "r") as inp, open(OUTPUT_FILE, "w") as out:
    current_chain = "A"
    last_residue = -1

    for line in inp:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            residue = int(line[22:26].strip())
            resname = line[17:20].strip()

            if resname == "SOL":
                break
            elif residue < last_residue:
                current_chain = increment_chain(current_chain)

            line = line[:21] + current_chain + line[22:]
            last_residue = residue
            out.write(line)

# Create output directory
if os.path.exists(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Split the output pdb file into separate files based on chain names
with open(OUTPUT_FILE, "r") as f:
    for line in f:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            chain = line[21]
            output_file = f"{OUTPUT_DIR}/chain_{chain}.pdb"
            with open(output_file, "a") as out_f:
                out_f.write(line)

def cleanup():
    files_to_remove = [
        "#*", "conf.gro", "editconf.pdb", "solvate.pdb", "ions.mdp", "ions.tpr", "index.ndx", "mdout.mdp",
        "topol.top", "posre.itp", 
        "complex*.pdb", "free*.pdb", "mol.pdb", "peptide_mutate.pdb",
        "topol_*","pmx_*", "posre_*"
    ]
    for pattern in files_to_remove:
        for file in glob.glob(pattern):
            Path(file).unlink(missing_ok=True)
    shutil.rmtree("chains")

# Remove previous topol.top file
Path("topol.top").unlink(missing_ok=True)

# Run gmx to get conf.gro
# subprocess.run(["gmx", "pdb2gmx", "-f", "chains/chain_C.pdb", "-ff", "charmm36-feb2021", "-water", "tip3p"], check=True)
subprocess.run(["gmx", "pdb2gmx", "-f", "chains/chain_C.pdb", "-ff", "charmm36m-mut", "-water", "tip3p"], check=True)

# Mutate conf.gro using pmx
subprocess.run(["pmx", "mutate", "-f", "conf.gro", "-ff", "charmm36m-mut", "-o", "peptide_mutate.pdb", "--script", "mut"])

# Divide (and fix the names of) complex PDB into multiple fragments using segname
tcl_script = '''
mol new peptide_mutate.pdb
set sel [atomselect top "all"]
$sel set chain X
$sel writepdb peptide_mutate.pdb
quit
'''

with open("tcl", "w") as tcl_file:
    tcl_file.write(tcl_script)

subprocess.run(["vmd", "-dispdev", "text", "-e", "tcl"])
Path("tcl").unlink(missing_ok=True)

# Get free.pdb
with open("peptide_mutate.pdb", "r") as f, open("mol.pdb", "w") as mol_f:
    for line in f:
        if line.startswith("ATOM"):
            mol_f.write(line)

subprocess.run(["gmx", "pdb2gmx", "-f", "mol.pdb", "-ff", "charmm36m-mut", "-water", "tip3p", "-o", "free.pdb"], check=True)

# Get complex.pdb (complex) and topol.top (complex)
for chain in ["A", "B"]:
    with open(f"chains/chain_{chain}.pdb", "r") as chain_f, open("mol.pdb", "a") as mol_f:
        for line in chain_f:
            if line.startswith("ATOM"):
                mol_f.write(line)
subprocess.run(["gmx", "pdb2gmx", "-f", "mol.pdb", "-ff", "charmm36m-mut", "-water", "tip3p", "-o", "complex.pdb"])

# Generate hybrid topology using pmx
subprocess.run(["pmx", "gentop", "-p", "topol.top", "-ff", "charmm36m-mut"])

with open("pmxtop.top", "r") as f, open("pmxtop_free.top", "w") as free_f:
    for line in f:
        if not (line.startswith("Protein_chain_A") or line.startswith("Protein_chain_B")):
            free_f.write(line)

# Create ions.mdp file
ions_mdp_content = '''
integrator  = md
dt          = 0.001
nsteps      = 100
nstlist         = 1
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
'''

with open("ions.mdp", "w") as ions_mdp_file:
    ions_mdp_file.write(ions_mdp_content)

# Remove and create directories
shutil.rmtree("free", ignore_errors=True)
shutil.rmtree("complex", ignore_errors=True)
Path("free/pdb2gmx").mkdir(parents=True, exist_ok=True)
Path("complex/pdb2gmx").mkdir(parents=True, exist_ok=True)

src_path = os.path.join(os.environ["GMXLIB"], "charmm36m-mut.ff")
os.symlink(src_path, os.path.join(f"free/pdb2gmx", "charmm36m-mut.ff"))
os.symlink(src_path, os.path.join(f"complex/pdb2gmx", "charmm36m-mut.ff"))

# Function for creating box, solvating, and adding ions
def prepare_system(system_type, input_pdb, output_top, box_size=None):
    # Create the box based on system_type
    if system_type == "complex":
        subprocess.run(f"echo 0 | gmx editconf -f {input_pdb} -o editconf.pdb -princ -d 1.0 -bt cubic", shell=True, check=True, text=True)
    elif system_type == "free":
        if not box_size:
            print("Box size must be provided for 'free'. Exiting.")
            sys.exit(1)
        subprocess.run(f"echo 0 | gmx editconf -f {input_pdb} -o editconf.pdb -box {box_size}", shell=True, check=True, text=True)

    # Solvate the box
    subprocess.run(["gmx", "solvate", "-cp", "editconf.pdb", "-o", "solvate.pdb", "-p", output_top])

    # Add ions to make the system neutral and of 0.15 M NaCl
    subprocess.run(["gmx", "grompp", "-f", "ions.mdp", "-c", "solvate.pdb", "-o", "ions.tpr", "-p", output_top, "-maxwarn", "1"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    subprocess.run(f"echo 13 | gmx genion -s ions.tpr -o {system_type}_ionized.pdb -conc 0.15 -neutral -p {output_top}", shell=True, check=True, text=True)

    # Create index file
    with subprocess.Popen(f"gmx make_ndx -f {system_type}_ionized.pdb", shell=True, stdin=subprocess.PIPE, text=True) as proc:
        proc.communicate(f"q\n")

    # Copy files to the respective directory
    files_to_copy = ["index.ndx"] + glob.glob("pmx_*") + glob.glob("posre_*")
    for file in files_to_copy:
        shutil.copy(file, f"{system_type}/pdb2gmx/{file}")

    shutil.move(output_top, f"{system_type}/pdb2gmx/topol.top")
    shutil.move(f"{system_type}_ionized.pdb", f"{system_type}/pdb2gmx/ionized.pdb")

    # Calculate the box size
    subprocess.run(["gmx", "editconf", "-f", f"{system_type}/pdb2gmx/ionized.pdb", "-o", "tmp.gro"])
    with open("tmp.gro", "r") as f:
        box_size = f.readlines()[-1]
    os.remove("tmp.gro")

    # Return the box size
    return box_size

# Prepare complex system
box_size_complex = prepare_system("complex", "complex.pdb", "pmxtop.top")

# Prepare free system
box_size_free = prepare_system("free", "free.pdb", "pmxtop_free.top", box_size_complex)

# Cleanup temporary files
cleanup()
