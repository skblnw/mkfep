import os
import sys
import subprocess
import math
import logging
from pathlib2 import Path
import argparse

# Setting up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Argument parsing
parser = argparse.ArgumentParser(description="Script description here")
parser.add_argument("--positions", type=int, nargs='+', default=[2, 3, 4, 5, 6, 7, 8],
                    help="Positions to process")
parser.add_argument("--residues", type=str, nargs='+', default=["ALA"],
                    help="Residues to process")
parser.add_argument("--exec_mbar", type=bool, default=True,
                    help="Execute MBAR or not")
parser.add_argument("--directories_num", type=int, default=16,
                    help="Number of directories")
args = parser.parse_args()

base = os.getcwd()
positions = args.positions
residues = args.residues
exec_mbar = args.exec_mbar
directories_num = args.directories_num

def check_command_exists(command):
    command_found = any(
        os.access(os.path.join(path, command), os.X_OK)
        for path in os.environ["PATH"].split(os.pathsep)
    )
    return command_found

def mkgmx_pymbar(base_dir):
    xvg_dir = base_dir / "xvg"
    xvg_dir.mkdir(parents=True, exist_ok=True)

    for i in range(directories_num):
        symlink_path = xvg_dir / 'dhdl.{}.xvg'.format(i)
        if symlink_path.exists() or symlink_path.is_symlink():
            os.remove(str(symlink_path))
        os.symlink(str(base_dir / 'dir{}/md.xvg'.format(i)), str(symlink_path))

    devnull = open(os.devnull, 'w')
    cmd = ['alchemical_analysis', '-d', str(xvg_dir), '-u', 'kcal', '-g', '-w', '-s', '1000']
    subprocess.run(cmd, stdout=devnull, stderr=subprocess.STDOUT)
    devnull.close()
    print("Done")

if check_command_exists("alchemical_analysis"):
    combined_resfile_path = "{}/combined_results.txt".format(base)
    if os.path.exists(combined_resfile_path):
        os.remove(combined_resfile_path)

    with open(combined_resfile_path, "a") as combined_resfile:
        for pos in positions:
            message = "Processing position: {} ".format(pos)
            print(message)
            for residue in residues:
                message = "  Processing residue: {} ".format(residue)
                print(message)
                combined_resfile.write("{} {}: ".format(pos, residue))
                
                free_dir = Path("{}/pos{}/{}/win{}.t1/free".format(base, pos, residue, directories_num))
                complex_dir = Path("{}/pos{}/{}/win{}.t1/complex".format(base, pos, residue, directories_num))
                
                if not free_dir.exists():
                    combined_resfile.write("Free state not exist; \n")
                    continue
                elif not complex_dir.exists():
                    combined_resfile.write("Bound state not exist; \n")
                    continue
                else:
                    if exec_mbar:
                        message = "    Free state ... "
                        sys.stdout.write(message)
                        sys.stdout.flush()
                        mkgmx_pymbar(free_dir)
                        message = "    Bound state ... "
                        sys.stdout.write(message)
                        sys.stdout.flush()
                        mkgmx_pymbar(complex_dir)
                    
                    free_result = free_dir / "xvg/results.txt"
                    complex_result = complex_dir / "xvg/results.txt"
                    
                    if not free_result.exists():
                        combined_resfile.write("Free state broken; \n")
                        continue
                    elif not complex_result.exists():
                        combined_resfile.write("Bound state broken; \n")
                        continue
                    else:
                        with open(str(complex_result), 'r') as f:
                            complex_total_line = next(line for line in f if "TOTAL:" in line)
                        with open(str(free_result), 'r') as f:
                            free_total_line = next(line for line in f if "TOTAL:" in line)
                        
                        b, bd = [float(value) for value in complex_total_line.split()[16:20:2]]
                        f, fd = [float(value) for value in free_total_line.split()[16:20:2]]
                        
                        result = b - f
                        err = math.sqrt(bd ** 2 + fd ** 2)
                        result_formatted = "{:7.2f}".format(result)
                        err_formatted = "{:4.2f}".format(err)
                        combined_resfile.write("{} +- {}\n".format(result_formatted, err_formatted))
else:
    logger.error("Command 'alchemical_analysis' not found. Please install it or check your system path.")
