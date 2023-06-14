import os
import sys
import subprocess
import shutil
import glob
from pathlib import Path

def BASH(cmd):
    subprocess.run(cmd, check=True, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

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

def check_required_files(required_files):
    for file in required_files:
        if not os.path.isfile(file):
            print(f"Error: Required file '{file}' not found.")
            sys.exit(1)

def convert_gro_to_pdb(gro_file, pdb_file):
    if os.path.exists(gro_file):
        subprocess.run(["gmx", "editconf", "-f", gro_file, "-o", pdb_file], check=True)
        return pdb_file
    elif os.path.exists(pdb_file):
        return pdb_file
    else:
        raise FileNotFoundError(f"Neither {pdb_file} nor {gro_file} exist.")

def assign_chain_identifiers(input_file, output_file):

    def increment_chain(chain):
        return chr(ord(chain) + 1)
        
    chain_data = {}
    with open(input_file, "r") as inp, open(output_file, "w") as out:
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

                # Update chain data
                if current_chain not in chain_data:
                    chain_data[current_chain] = {'first_residue': residue, 'last_residue': residue}
                else:
                    chain_data[current_chain]['last_residue'] = residue

    # Log chain data
    print(f"Assigned {len(chain_data)} chains.")
    for chain, data in chain_data.items():
        print(f"Chain {chain}: first residue = {data['first_residue']}, last residue = {data['last_residue']}")

def create_output_directory(output_dir):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

def split_pdb_file(output_file, output_dir):
    with open(output_file, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21]
                output_file = f"{output_dir}/chain_{chain}.pdb"
                with open(output_file, "a") as out_f:
                    out_f.write(line)

# Variables
PREFIX = "md"
OUTPUT_FILE = f"{PREFIX}_with_chains.pdb"
OUTPUT_DIR = "chains"
required_files = ["mut"]

# Check if required files exist
check_required_files(required_files)

# Convert gro to pdb if needed
INPUT_FILE = convert_gro_to_pdb(f"{PREFIX}.gro", f"{PREFIX}.pdb")

# Assign chain identifiers to repeating residue number ranges
assign_chain_identifiers(INPUT_FILE, OUTPUT_FILE)

# Create output directory
create_output_directory(OUTPUT_DIR)

# Split the output pdb file into separate files based on chain names
split_pdb_file(OUTPUT_FILE, OUTPUT_DIR)

def generate_topology():
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

def create_directory(path):
    os.makedirs(path, exist_ok=True)

def write_to_file(path, content):
    with open(path, 'w') as f:
        f.write(content)

def copy_file(src, dest):
    shutil.copyfile(src, dest)

def change_directory(path):
    os.chdir(path)

def NTERIO(gmxtopfn, pmxtopfn, outfn, m=4):
  with open(gmxtopfn, 'r') as fg:
    og=[[]]
    n = 0
    for l in fg.readlines():
      if l.startswith('['):
        n=n+1
        og.append([])
      if n in [2]:
        if any([l.startswith('    %d'%(i+1),1)  for i in range(6)]) : og[n].append(l)
      if n in [3,4]:
        if any([l.startswith('    %d'%(i+1),p)  for p in [0,6]  for i in range(m)]) : og[n].append(l)
      if n in [5]:
        if any([l.startswith('    %d'%(i+1),p)  for p in [0,6,12]  for i in range(m)]) : og[n].append(l)
      if n in [6,7]:
        if any([l.startswith('    %d'%(i+1),p)  for p in [0,6,12,18]  for i in range(m)]) : og[n].append(l)
      if n in [8]:
        if any([l.startswith('    %d'%(i+1),p)  for p in [0,6,12,18,24]  for i in range(m)]) : og[n].append(l)
  with open(pmxtopfn, 'r') as fp:
    out = []
    n = 0
    for l in fp.readlines():
      if l.startswith('['):
        n=n+1
      if n<2 or n>8:
        out.append(l)
      if n in [2]:
        if not any([l.startswith('    %d'%(i+1),1)  for i in range(6)]) : out.append(l)
      if n in [3,4]:
        if not any([l.startswith('     %d'%(i+1),p)  for p in [0,7]  for i in range(m)]) : out.append(l)
      if n in [5]:
        if not any([l.startswith('     %d'%(i+1),p)  for p in [0,7,14]  for i in range(m)]) : out.append(l)
      if n in [6,7]:
        if not any([l.startswith('     %d'%(i+1),p)  for p in [0,7,14,21]  for i in range(m)]) : out.append(l)
      if n in [8]:
        if not any([l.startswith('     %d'%(i+1),p)  for p in [0,7,14,21,28]  for i in range(m)]) : out.append(l)
      if l.startswith('['):
        out = out + og[n]
  with open(outfn,'w') as f:
    f.writelines(out)

def generate_topology_new(nter=False):
    # Remove previous topol.top file
    Path("topol.top").unlink(missing_ok=True)

    # Run gmx to get conf.gro
    BASH('gmx pdb2gmx -f chains/chain_C.pdb -o original.pdb -ff charmm36m-mut -water tip3p')

    if nter:
        BASH("grep ' [34]  H[23]  ' original.pdb > znt.pdb")
        BASH("sed 's/ 2  H1  / 2  HN  /' original.pdb -i")
        BASH("sed 's/H2.*C   1/C   ZNT C   0/' znt.pdb -i")
        BASH("sed 's/H3.*C   1/O   ZNT C   0/' znt.pdb -i")

    # Mutate conf.gro using pmx
    BASH('pmx mutate -f original.pdb -o peptide_mutate.pdb -ff charmm36m-mut --script mut')

    BASH("sed '/^[^A]/d' peptide_mutate.pdb -i")

    BASH('gmx pdb2gmx -f peptide_mutate.pdb -o free.pdb -ff charmm36m-mut -water tip3p')

    BASH('cat peptide_mutate.pdb chains/chain_A.pdb chains/chain_B.pdb > mutant.pdb')

    BASH('gmx pdb2gmx -f mutant.pdb -o complex.pdb -ff charmm36m-mut -water tip3p')

    if nter:
        # Copy files
        for file in ['topol_Protein_chain_C.itp', 'posre_Protein_chain_C.itp']:
            copy_file(file, f'{file}.copy')

        BASH('cat znt.pdb peptide_mutate.pdb chains/chain_A.pdb chains/chain_B.pdb > mutant.pdb')

        BASH('echo 3 0 0 0 0 0 | gmx pdb2gmx -f mutant.pdb -o conf.pdb -ff charmm36m-mut -water tip3p -ter')

        BASH('pmx gentop -p topol.top -ff charmm36m-mut')

        # Copy pmx_topol_Protein_chain_C.itp
        copy_file('pmx_topol_Protein_chain_C.itp', 'pmx_topol_Protein_chain_C.itp.copy')

        # Run NTERIO
        NTERIO('topol_Protein_chain_C.itp.copy', 'pmx_topol_Protein_chain_C.itp.copy', 'pmx_topol_Protein_chain_C.itp')

        for file in ['topol_Protein_chain_C.itp', 'posre_Protein_chain_C.itp']:
            copy_file(f'{file}.copy', file)
    else:
        # Generate hybrid topology using pmx
        BASH('pmx gentop -p topol.top -ff charmm36m-mut')

    # Write to pmxtop_free.top
    with open("pmxtop.top", "r") as f, open("pmxtop_free.top", "w") as free_f:
        for line in f:
            if not (line.startswith("Protein_chain_A") or line.startswith("Protein_chain_B")):
                free_f.write(line)

def check_terminal_mutations(file):
    nter = False
    cter = False
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('1'):
                nter = True
            elif line.startswith('9'):
                cter = True
    return nter, cter

nter, cter = check_terminal_mutations('mut')

generate_topology_new(nter=nter)

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

def create_directory(path):
    os.makedirs(path, exist_ok=True)

def write_to_file(path, content):
    with open(path, 'w') as f:
        f.write(content)

def copy_file(src, dest):
    shutil.copyfile(src, dest)

def change_directory(path):
    os.chdir(path)

def BASH(cmd):
    subprocess.run(cmd, check=True, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def prepare_system_new(system_type, input_pdb, output_top, box_size=None):
    # Create the box based on system_type
    if system_type == "complex":
        BASH("echo 0 | gmx editconf -f %s -o editconf.pdb -princ -d 1.0 -bt cubic" % input_pdb)
    elif system_type == "free":
        if not box_size:
            print("Box size must be provided for 'free'. Exiting.")
            sys.exit(1)
        BASH("echo 0 | gmx editconf -f %s -o editconf.pdb -box %s" % (input_pdb, box_size))

    # Solvate the box
    BASH("gmx solvate -cp editconf.pdb -o solvate.pdb -p %s" % output_top)

    # Add ions to make the system neutral and of 0.15 M NaCl
    BASH("gmx grompp -f ions.mdp -c solvate.pdb -o ions.tpr -p %s -maxwarn 1" % output_top)
    BASH("echo 13 | gmx genion -s ions.tpr -o %s_ionized.pdb -conc 0.15 -neutral -p %s" % (system_type, output_top))

    # Create index file
    with subprocess.Popen("gmx make_ndx -f %s_ionized.pdb" % system_type, shell=True, stdin=subprocess.PIPE, text=True) as proc:
        proc.communicate(f"q\n")

    # Copy files to the respective directory
    files_to_copy = ["index.ndx"] + glob.glob("pmx_*") + glob.glob("posre_*")
    for file in files_to_copy:
        shutil.copy(file, "%s/pdb2gmx/%s" % (system_type, file))

    shutil.move(output_top, "%s/pdb2gmx/topol.top" % system_type)
    shutil.move("%s_ionized.pdb" % system_type, "%s/pdb2gmx/ionized.pdb" % system_type)

    # Calculate the box size
    BASH("gmx editconf -f %s/pdb2gmx/ionized.pdb -o tmp.gro" % system_type)
    with open("tmp.gro", "r") as f:
        box_size = f.readlines()[-1]
    os.remove("tmp.gro")

    # Return the box size
    return box_size

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
    subprocess.run(f"echo SOL | gmx genion -s ions.tpr -o {system_type}_ionized.pdb -conc 0.15 -neutral -p {output_top}", shell=True, check=True, text=True)

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
