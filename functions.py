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



def copy_atoms(src, dest):
    with open(src, 'r') as src_f, open(dest, 'w') as dest_f:
        for line in src_f:
            if line.startswith("ATOM"):
                dest_f.write(line)




def generate_topology_zhong():
    # Define paths
    free_path = f'{args.mutlabel}/free/pdb2gmx'
    complex_path = f'{args.mutlabel}/complex/pdb2gmx'
    pmx_script_path = f'{args.mutlabel}/pmx_script'

    # Create directories
    create_directory(free_path)
    create_directory(complex_path)

    # Write to pmx script
    write_to_file(pmx_script_path, tmppmxs)

    # Check for terminal mutations
    nter, cter = '1' in tmppmxs, f'{len(args.mutseq)}' in tmppmxs

    # Change to free path
    change_directory(f'{rpath}/{free_path}')

    # Copy conf.C.pdb to original.pdb
    copy_file(f'{rpath}/conf.C.pdb', 'original.pdb')

    # Handle N-terminal mutation
    if nter:
        BASH("grep 'H[23]' original.pdb > znt.pdb")
        BASH("sed '/H1/s/H1/HN/' original.pdb -i")
        BASH("sed 's/H2.*C   1/C   ZNT C   0/' znt.pdb -i")
        BASH("sed 's/H3.*C   1/O   ZNT C   0/' znt.pdb -i")

    # Run pmx mutate
    BASH(f'pmx mutate -f original.pdb -ff {args.ff} -o mutant.pdb --script {rpath}/{pmx_script_path}')

    # Clean up mutant.pdb
    BASH("sed -i '/^[^A]/d' mutant.pdb")

    # Change to complex path
    change_directory(f'{rpath}/{complex_path}')

    # Generate complex mutant structure
    BASH(f'cat {rpath}/conf.A.pdb {rpath}/conf.B.pdb {rpath}/{free_path}/mutant.pdb > mutant.pdb')

    # Run gmx pdb2gmx
    BASH(f'gmx pdb2gmx -f mutant.pdb -ff {args.ff} -water tip3p -o conf.pdb')

    # Handle N-terminal mutation
    if nter:
        # Copy files
        for file in ['conf.pdb', 'topol_Protein_chain_C.itp', 'posre_Protein_chain_C.itp']:
            copy_file(file, f'{file}.copy')

        # Generate mutant.pdb
        BASH(f'cat {rpath}/conf.A.pdb {rpath}/conf.B.pdb {rpath}/{free_path}/znt.pdb {rpath}/{free_path}/mutant.pdb > mutant.pdb')

        # Run gmx pdb2gmx
        BASH(f'echo 0 0 0 0 3 0 | gmx pdb2gmx -f mutant.pdb -ff {args.ff} -water tip3p -o conf.pdb -ter')

        # Run pmx gentop
        BASH(f'pmx gentop -p topol.top -ff {args.ff}')

        # Copy pmx_topol_Protein_chain_C.itp
        copy_file('pmx_topol_Protein_chain_C.itp', 'pmx_topol_Protein_chain_C.itp.copy')

        # Run NTERIO
        NTERIO('topol_Protein_chain_C.itp.copy', 'pmx_topol_Protein_chain_C.itp.copy', 'pmx_topol_Protein_chain_C.itp')

        # Move files
        for file in ['topol_Protein_chain_C.itp', 'posre_Protein_chain_C.itp', 'conf.pdb']:
            os.rename(f'{file}.copy', file)
    else:
        # Run pmx gentop
        BASH(f'pmx gentop -p topol.top -ff {args.ff}')







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

    # Divide (and fix the names of) complex PDB into multiple fragments using segname
    # tcl_script = '''
    # mol new peptide_mutate.pdb
    # set sel [atomselect top "all"]
    # $sel set chain X
    # $sel writepdb peptide_mutate.pdb
    # quit
    # '''
    # write_to_file("tcl", tcl_script)
    # BASH('vmd -dispdev text -e tcl')
    # Path("tcl").unlink(missing_ok=True)

    # # Get free.pdb
    # copy_atoms("peptide_mutate.pdb", "mol.pdb")
    # BASH('gmx pdb2gmx -f mol.pdb -ff charmm36m-mut -water tip3p -o free.pdb')

    # Get complex.pdb (complex) and topol.top (complex)
    # for chain in ["A", "B"]:
    #     copy_atoms(f"chains/chain_{chain}.pdb", "mol.pdb")

    BASH("sed '/^[^A]/d' peptide_mutate.pdb -i")

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


