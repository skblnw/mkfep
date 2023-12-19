"""
A simple tool to automatically get results from NAMD simulation output. Author: Xufan Gao, Kevin C. Chan

Test: at `Lei2Guanqiao`

"""
__version__ = 0.0

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import os, glob, argparse, logging
import numpy as np
import pandas as pd
import alchemlyb
from alchemlyb.parsing import namd
from alchemlyb.estimators import MBAR, BAR, TI
methods = {'BAR': BAR, 'MBAR': MBAR, 'TI': TI}
from alchemlyb.postprocessors.units import to_kcalmol, R_kJmol, kJ2kcal
k_b = R_kJmol * kJ2kcal
from alchemlyb.visualisation import plot_convergence
from alchemlyb.convergence import forward_backward_convergence


# user IO
class CustomMetavarFormatter(argparse.RawTextHelpFormatter):
    """
    Reference: https://devpress.csdn.net/python/62fe2a1dc67703293080479b.html

    If the optional takes a value, format is: ``-s ARGS, --long ARGS``; Now changed to ``-s, --long ARGS``
    """

    def _format_action_invocation(self, action):
        if not action.option_strings:
            metavar, = self._metavar_formatter(action, action.dest)(1)
            return metavar
        else:
            parts = []
            if action.nargs == 0:
                parts.extend(action.option_strings)
            else:
                default = action.dest.upper()
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    # parts.append('%s %s' % (option_string, args_string))
                    parts.append('%s'%option_string)
                parts[-1] += ' %s'%args_string
            return ', '.join(parts)


def define_parser(version):
    usage = "python fepout-analysis.py [-h] -f xxx.csv [options]"
    des = "Automated BAR analysis for NAMD fepout files, considering both forward and backward runs.\n\n"
    des += "An example input csv file should look like (empty lines and lines starting with # are commented out):\n"
    des += "forward,bound,1.csv\nforward,bound,2.csv\n...\nbackward,free,3.csv\nbackward,free,4.csv\n"
    epilog = 'Welcome to star this project at https://github.com/gxf1212/FEbuilder!'

    parser = argparse.ArgumentParser(prog='fepout-analysis', description=des, usage=usage, epilog=epilog,
                                     formatter_class=CustomMetavarFormatter)
    basic = parser.add_argument_group('Basic')
    basic.add_argument('-f', '--fepout', type=str, metavar='FILE', required=True,
                       help='A csv file containing paths to the fepout files. Required.')
    basic.add_argument('-T', '--temperature', type=float, metavar='TEMP', default=310,
                        help='The temperature (in Kelvin) of the simulation. Default is 310.')
    basic.add_argument('-m', '--method', metavar='ESTIMATOR', choices=methods.keys(), default='BAR',
                        help='Estimator to calculate free energy differences. Default is BAR.')
    basic.add_argument('-l', '--logfile', type=str, metavar='LOG', default='fepout-analysis.log',
                       help='Path to the log file. Default is fepout-analysis.log.')
    # advanced = parser.add_argument_group('Advanced')
    # not sure yet, basically options for plotting
    return parser


def parse_input(parser):
    args = parser.parse_args()
    args.method = methods[args.method]
    args.logger = define_logger(args.logfile)
    return args


def define_logger(logfile):
    # copied from FEbuilder
    # Create a logger object
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    # Create a handler for writing log messages to a file
    logfile = os.path.join(os.path.abspath('.'), logfile)
    try:
        os.remove(logfile)  # must remove the old before creation of logger in this run
    except FileNotFoundError as e:
        pass
    file_handler = logging.FileHandler(logfile)
    file_handler.setLevel(logging.INFO)
    # Create a handler for writing log messages to the console
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    # Add the handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    return logger


# Preprocessing, and tools
def namd_preprocess_fepout(u_nk_fwd, u_nk_rev, edit_lambda=True):
    """
    Sort by lambda value and timesteps; unify column names, especially transform backward lambdas to 1 - lambda
    Avoid numerical precision errors for very small numbers. Nobody will set lambda to 1e-8....

    Args: both forward and backward DataFrames

    Returns: revised DataFrames
    """
    u_nk_fwd = u_nk_fwd.sort_index(level=u_nk_fwd.index.names[1:])
    u_nk_fwd.columns = [round(c, 8) for c in u_nk_fwd.columns]
    # reverse lambda values
    if edit_lambda:
        u_nk_rev.columns = [round(1-c, 8) for c in u_nk_rev.columns]  # Update column names
        u_nk_rev.index = u_nk_rev.index.set_levels([round(1-c, 9) for c in u_nk_rev.index.levels[1]], level=1)  # Update fep-lambda indices
        u_nk_rev = u_nk_rev.sort_index(level=u_nk_fwd.index.names[1:]).sort_index(axis=1)  # sort
    return u_nk_fwd, u_nk_rev


def namd_fit_fwbk(u_nk_fwd, u_nk_rev, method=None):
    """
    Combine values in forward and backward DataFrames. For column lambda_i, collect values for lambda_i-1 and lambda_i+1 (except 0 and 1).

    Return the combined DataFrame or the result.
    """
    # combine the two
    u_nk_fwd.replace(0, np.nan, inplace=True)
    u_nk_fwd[u_nk_fwd.isnull()] = u_nk_rev
    u_nk_fwd.replace(np.nan, 0, inplace=True)
    u_nk = u_nk_fwd.sort_index(level=u_nk_fwd.index.names[1:])
    # sort final dataframe by `fep-lambda` (as opposed to `timestep`)
    if method is None:  # just combine forward/backward
        return u_nk
    else:  # fit
        bar = method().fit(u_nk)
        delta_f = to_kcalmol(bar.delta_f_)
        return delta_f


def calculate_sub_diagonal_sqrt_sum(df):
    """
    Return the stderr of each window, and the square root of the sum of squares of these elements.
    """
    # Ensure the DataFrame is square
    if df.shape[0] != df.shape[1]:
        raise ValueError("DataFrame must be square")
    # Calculate the sum of squares of the sub-diagonal elements
    sub_diagonal_elements = np.array([df.iloc[i+1, i] for i in range(df.shape[0]-1)])
    root_sum_squares = np.linalg.norm(sub_diagonal_elements)
    return sub_diagonal_elements, root_sum_squares


def extract_vdw_elec_dataframes(filename, T):
    """
    Extracts vdW and Elec values from a given FEP output file and returns two DataFrames.

    Args:
    filename (str): The path to the FEP output file.

    Returns:
    tuple of DataFrame: Two pandas DataFrames containing vdW and Elec values. Null values filled with 0.
    """
    dict_u_nk_vdw = {}
    dict_u_nk_elec = {}
    lambda_values = set()
    time_lambda_pairs = set()
    parsing = False
    beta = 1 / (k_b * T)

    with open(filename, 'r') as file:
        for line in file:
            if "#Free energy change for lambda window" in line:
                parsing = False
            if parsing and line.startswith("FepEnergy:"):
                # collect data
                parts = line.split()
                time = float(parts[1])
                elec_l, elec_ldl = float(parts[2]), float(parts[3])
                vdw_l, vdw_ldl = float(parts[4]), float(parts[5])
                elec_diff = round(elec_ldl - elec_l, 8)
                vdw_diff = round(vdw_ldl - vdw_l, 8)
                # collect indices
                dict_u_nk_elec[(time, lambda1)] = dict_u_nk_elec.get((time, lambda1), {})
                dict_u_nk_vdw[(time, lambda1)] = dict_u_nk_vdw.get((time, lambda1), {})
                dict_u_nk_elec[(time, lambda1)][lambda2] = elec_diff
                dict_u_nk_vdw[(time, lambda1)][lambda2] = vdw_diff
                time_lambda_pairs.add((time, lambda1))
                lambda_values.add(lambda1)
                lambda_values.add(lambda2)
            if line.startswith("#NEW FEP WINDOW:"):
                # set the current lambda values
                parts = line.split()
                lambda1, lambda2 = round(float(parts[-3]), 8), round(float(parts[-1]), 8)
            if "STEPS OF EQUILIBRATION AT LAMBDA" in line:
                parsing = True

    # Add the last lambda value to dictionary, and fill null with 0
    for time in set(time for time, _ in time_lambda_pairs):
        time_lambda_pairs.add((time, lambda2))
        dict_u_nk_elec.setdefault((time, lambda2), {}).update({l: 0 for l in lambda_values})
        dict_u_nk_vdw.setdefault((time, lambda2), {}).update({l: 0 for l in lambda_values})

    # Create DataFrame structure
    sorted_lambdas = sorted(list(lambda_values))
    index = pd.MultiIndex.from_tuples(time_lambda_pairs, names=['time', 'fep-lambda'])
    u_nk_elec = pd.DataFrame(index=index, columns=sorted_lambdas).fillna(0)
    u_nk_vdw = pd.DataFrame(index=index, columns=sorted_lambdas).fillna(0)
    u_nk_elec.attrs = u_nk_vdw.attrs = {'temperature': T, 'energy_unit': 'kT'}  # otherwise cannot be processed by pymbar

    # Populate DataFrames
    for (time, lambda1), values in dict_u_nk_elec.items():
        for lambda2, elec_value in values.items():
            u_nk_elec.at[(time, lambda1), lambda2] = float(elec_value) * beta
    for (time, lambda1), values in dict_u_nk_vdw.items():
        for lambda2, vdw_value in values.items():
            u_nk_vdw.at[(time, lambda1), lambda2] = float(vdw_value) * beta

    return u_nk_elec.sort_index(level=u_nk_elec.index.names[1:]), u_nk_vdw.sort_index(level=u_nk_vdw.index.names[1:])


# fitting for free energy change
def update_results(result, stderr, prefix, fpath, bpath, method, T):
    """
    Sub-process for get_results.

    Parameters:
    result (dict): Dictionary to store results.
    stderr (dict): Dictionary to store standard errors.
    prefix (str): Prefix to label the results.
    fpath (list): List of file paths for forward calculations.
    bpath (list): List of file paths for backward calculations.
    method (callable): Method for free energy calculation.
    T (float): Temperature for the calculation.
    """

    def process_dict(f_df, b_df, component_label):
        u_nk_fwd, u_nk_rev = namd_preprocess_fepout(
            alchemlyb.concat(f_df),
            alchemlyb.concat(b_df)
        )
        u_nk = namd_fit_fwbk(u_nk_fwd, u_nk_rev, method=None)
        bar = method().fit(u_nk)
        delta_f = to_kcalmol(bar.delta_f_)
        ddeltaf = to_kcalmol(bar.d_delta_f_)

        result[prefix][component_label] = delta_f[1.0][0.0]
        stderr[prefix][component_label] = calculate_sub_diagonal_sqrt_sum(ddeltaf)[1]

    # Total
    f_df = [namd.extract_u_nk(fp, T) for fp in fpath]
    b_df = [namd.extract_u_nk(bp, T) for bp in bpath]
    process_dict(f_df, b_df, 'dG')
    args.logger.info(f'Done {prefix} dG')

    # Electrostatics and van der Waals components
    df_fwd = [extract_vdw_elec_dataframes(fp, T) for fp in fpath]
    df_rev = [extract_vdw_elec_dataframes(bp, T) for bp in bpath]
    for component, idx in [('elec', 0), ('vdw', 1)]:
        process_dict(
            [d[idx] for d in df_fwd],
            [d[idx] for d in df_rev],
            component
        )
    args.logger.info(f'Done {prefix} decomposition')

    return result, stderr


def get_results(filelist):
    """
    Get the results from the fepout files, handling multiple file paths for each key.

    Parameters:
    filelist (dict): Dictionary containing nested dictionaries for 'bound' and 'free' states, 
                     each with lists of file paths for 'forward' and 'backward'.

    Returns:
    tuple: Two dictionaries containing results and standard errors for bound and free states.
    """

    result = {}
    stderr = {}
    # Iterate over states ('bound' and 'free')
    for state in filelist:
        result[state] = {}
        stderr[state] = {}
        if 'forward' in filelist[state] and 'backward' in filelist[state]:
            fpath = filelist[state]['forward']  # list of file paths for forward calculations
            bpath = filelist[state]['backward']  # list of file paths for backward calculations
            update_results(result, stderr, state, fpath, bpath, args.method, args.temperature)
    return result, stderr


def read_file_list():
    """
    Reads and formats a combined CSV file for both bound and free states.
    Processes only lines with 'forward' or 'backward' keys and 'bound' or 'free' states.
    Skips lines that do not fit these criteria or start with '#'.

    Parameters:
    filepath (str): Path to the combined CSV file.

    Returns:
    dict: A dictionary containing 'bound' and 'free' states, 
          each with 'forward' and 'backward' keys mapped to lists of file paths.
    """

    filelist = {'bound': {'forward': [], 'backward': []}, 'free': {'forward': [], 'backward': []}}
    expected_states = {'bound', 'free'}
    expected_keys = {'forward', 'backward'}

    with open(args.fepout, 'r') as file:
        for line in file:
            # Skip comments and empty lines
            if line.startswith('#') or not line.strip():
                continue
            state, key, path = line.strip().split(',')
            if state in expected_states and key in expected_keys:
                filelist[state][key].append(path)
            else:
                # Optionally, you can print a warning or log skipped lines
                args.logger.info(f"Skipping line: {line.strip()}")

    return filelist


def output_results(result, stderr):
    """
    Print the results
    """
    args.logger.info('\nResult:')
    args.logger.info("Forward ddG is {:6.3f} ± {:.3f} kcal/mol".format(result['bound']['dG']-result['free']['dG'], np.linalg.norm([stderr['bound']['dG'], stderr['free']['dG']])))
    args.logger.info("ddG(elec)   is {:6.3f} ± {:.3f} kcal/mol".format(result['bound']['elec']-result['free']['elec'], np.linalg.norm([stderr['bound']['elec'], stderr['free']['elec']])))
    args.logger.info("ddG(vdw)    is {:6.3f} ± {:.3f} kcal/mol".format(result['bound']['vdw']-result['free']['vdw'], np.linalg.norm([stderr['bound']['vdw'], stderr['free']['vdw']])))


# wrap-up
def main():
    filelist = read_file_list()
    result, stderr = get_results(filelist)
    output_results(result, stderr)


global args  # more accessible to arguments
if __name__ == "__main__":
    parser = define_parser(__version__)
    args = parse_input(parser)
    main()


# debug
# if __name__ == "__main__":
#     path = '/data/work/make_hybrid_top/tests/Lei2Guanqiao'
#     os.chdir(path)
#     data = {
#         'fepout': os.path.join(path, 'D3N-N3D.csv'),
#         'method': BAR,
#         'temperature': 310
#         'logfile': 'fepout-analysis.log',
#     }
#     args = argparse.Namespace(**data)
#     args.logger = define_logger()
#     main()

