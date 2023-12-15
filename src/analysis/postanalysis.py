"""
A simple tool to automatically get results from NAMD simulation output.

"""


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import numpy as np
import os, glob
import pandas as pd
import alchemlyb
from alchemlyb.parsing import namd
from alchemlyb.estimators import MBAR, BAR
from alchemlyb.visualisation import plot_convergence
from alchemlyb.convergence import forward_backward_convergence
from alchemlyb.postprocessors.units import to_kcalmol, R_kJmol, kJ2kcal
from alchemlyb.visualisation import plot_mbar_overlap_matrix
k_b = R_kJmol * kJ2kcal


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
    # Ensure the DataFrame is square
    if df.shape[0] != df.shape[1]:
        raise ValueError("DataFrame must be square")
    # Calculate the sum of squares of the sub-diagonal elements
    sum_of_squares = sum(df.iloc[i+1, i]**2 for i in range(df.shape[0] - 1))
    return np.sqrt(sum_of_squares)  # Return the square root of the sum of squares


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


path = '/data/work/make_hybrid_top/tests/Lei2Guanqiao'
fwd = os.path.join(path, 'D3N')
bwd = os.path.join(path, 'N3D')
T = 310
method = BAR
max_rep = 6
result = {}
stderr = {}

for prefix in ['bound', 'free']:
    fpath = [f'{fwd}/{prefix}/t{repeat}/alchemy.fepout' for repeat in range(1, max_rep+1)]
    bpath = [f'{bwd}/{prefix}/t{repeat}/alchemy.fepout' for repeat in range(1, max_rep+1)]
    # total dG
    u_nk_fwd, u_nk_rev = namd_preprocess_fepout(alchemlyb.concat([namd.extract_u_nk(fp, T) for fp in fpath]),
                                                alchemlyb.concat([namd.extract_u_nk(bp, T) for bp in bpath]))
    u_nk = namd_fit_fwbk(u_nk_fwd, u_nk_rev, method=None)
    bar = method().fit(u_nk)
    delta_f = to_kcalmol(bar.delta_f_)
    result[prefix] = {'dG': delta_f[1.0][0.0]}
    ddeltaf = to_kcalmol(bar.d_delta_f_)
    stderr[prefix] = {'dG': calculate_sub_diagonal_sqrt_sum(ddeltaf)}
    print(f'Done {prefix} dG')
    # retrive elec/vdw from fepout
    df_fwd = [extract_vdw_elec_dataframes(fp, T) for fp in fpath]
    df_rev = [extract_vdw_elec_dataframes(bp, T) for bp in bpath]
    print(f'Done {prefix} decomp data')
    # elec
    u_nk_fwd_elec, u_nk_rev_elec = namd_preprocess_fepout(alchemlyb.concat([d[0] for d in df_fwd]),
                                                          alchemlyb.concat([d[0] for d in df_rev]))
    u_nk_elec = namd_fit_fwbk(u_nk_fwd_elec, u_nk_rev_elec, method=None)
    bar_ele = method().fit(u_nk_elec)
    delta_f = to_kcalmol(bar_ele.delta_f_)
    result[prefix]['elec'] = delta_f[1.0][0.0]
    ddeltaf = to_kcalmol(bar_ele.d_delta_f_)
    stderr[prefix]['elec'] = calculate_sub_diagonal_sqrt_sum(ddeltaf)
    print(f'Done {prefix} elec')
    # vdw
    u_nk_fwd_vdw, u_nk_rev_vdw = namd_preprocess_fepout(alchemlyb.concat([d[1] for d in df_fwd]),
                                                        alchemlyb.concat([d[1] for d in df_rev]))
    u_nk_vdw = namd_fit_fwbk(u_nk_fwd_vdw, u_nk_rev_vdw, method=None)
    bar_vdw = method().fit(u_nk_vdw)
    delta_f = to_kcalmol(bar_vdw.delta_f_)
    result[prefix]['vdw'] = delta_f[1.0][0.0]
    ddeltaf = to_kcalmol(bar_vdw.d_delta_f_)
    stderr[prefix]['vdw'] = calculate_sub_diagonal_sqrt_sum(ddeltaf)
    print(f'Done {prefix} vdw')
    # result[prefix]['residual'] = result[prefix]['dG'] - result[prefix]['elec'] - result[prefix]['vdw']

print()
print("Forward ddG is {:6.3f} ± {:.3f} kcal/mol".format(result['bound']['dG']-result['free']['dG'], np.linalg.norm([stderr['bound']['dG'], stderr['free']['dG']])))
print("ddG(elec)   is {:6.3f} ± {:.3f} kcal/mol".format(result['bound']['elec']-result['free']['elec'], np.linalg.norm([stderr['bound']['elec'], stderr['free']['elec']])))
print("ddG(vdw)    is {:6.3f} ± {:.3f} kcal/mol".format(result['bound']['vdw']-result['free']['vdw'], np.linalg.norm([stderr['bound']['vdw'], stderr['free']['vdw']])))
