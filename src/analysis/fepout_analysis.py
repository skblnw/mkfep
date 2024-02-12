"""
A simple tool to automatically get results from NAMD simulation output. Author: Xufan Gao; Acknowledgement: Kevin C. Chan

Test: at `Lei2Guanqiao`: python fepout_analysis.py -f D3N-N3D.csv -r D3N-N3D -c time-conv -dt 1.0
"""
__version__ = 0.0
expected_states = ['bound', 'free']
expected_keys = ['forward', 'backward']
expected_components = ['elec', 'vdw']

import warnings
warnings.simplefilter(action='ignore', category=Warning)
import os, shutil, sys, argparse, logging, collections
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties as FP
import alchemlyb
import alchemlyb_parsing_namd as apn
from alchemlyb.estimators import MBAR, BAR, TI
methods = {'BAR': BAR, 'MBAR': MBAR, 'TI': TI}
from alchemlyb.postprocessors.units import to_kcalmol, R_kJmol, kJ2kcal, get_unit_converter
k_b = R_kJmol*kJ2kcal
# from alchemlyb.visualisation import plot_convergence
# from alchemlyb.convergence import forward_backward_convergence


## user IO
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


def parse_arguments(version):
    usage = "python fepout_analysis.py [-h] -f xxx.csv [options]"
    des = "Automated BAR analysis for NAMD fepout files, considering both forward and backward runs.\n\n"
    des += "An example input csv file should look like (empty lines and lines starting with # are commented out):\n"
    des += "forward,bound,1.csv\nforward,bound,2.csv\n...\nbackward,free,3.csv\nbackward,free,4.csv\n"
    epilog = 'Welcome to star this project at https://github.com/gxf1212/FEbuilder!'

    parser = argparse.ArgumentParser(prog='fepout-analysis', description=des, usage=usage, epilog=epilog,
                                     formatter_class=CustomMetavarFormatter)
    parser.add_argument('-V', '--version', action='version', version='%(prog)s version '+str(version))

    basic = parser.add_argument_group('Basic')
    basic.add_argument('-f', '--fepout', type=str, metavar='FILE', required=True,
                       help='A csv file containing paths to the fepout files. Required.')
    basic.add_argument('-T', '--temperature', type=float, metavar='TEMP', default=310,
                       help='The temperature (in Kelvin) of the simulation. Default is 310.')
    basic.add_argument('-m', '--method', metavar='ESTIMATOR', choices=methods.keys(), default='BAR',
                       help='Estimator to calculate free energy differences. Default is BAR.')
    result = parser.add_argument_group('Result')
    result.add_argument('-r', '--result', type=str, metavar='DIR', default='fepout-result',
                        help='Result directory. Default is "fepout-result".')
    result.add_argument('-l', '--logfile', type=str, metavar='LOG', default='fepout-analysis.log',
                        help='Path to the log file. Default is fepout-analysis.log.')
    result.add_argument('-fo', '--fodgfile', type=str, metavar='dG', default='dg-fepout',
                        help='Pefix of the dG-from-fepout summary. Default is dg-fepout.')
    result.add_argument('-nd', '--no-decompose', default=False, action="store_true",
                         help="Do not do decomposition. Default is False.")
    result.add_argument('-d', '--dglfile', type=str, metavar='dG-L', default='dg-lambda',
                        help='Prefix of the dG-lambda results. Default is dg-lambda.')
    result.add_argument('-c', '--convfile', type=str, metavar='COV', default='',
                        help='Prefix of the time-convergence results, e.g. "conv".\nDefault is an empty string, i.e. analysis won\'t be done.')
    result.add_argument('-cn', '--conv-num', type=int, metavar='NUM', default=10,
                        help='Number of points in time-convergence analysis. Default is 10')
    result.add_argument('-dt', type=float, default=2.0, metavar='\u0394t',
                        help='NAMD parameter "timestep". Equilsteps or Numofsteps * deltaT = simulation time.\n'+
                             'Default is 2.0, unit is femtosecond (fs). For time-convergence analysis.')
    result.add_argument('-bs', '--bsfile', type=str, metavar='BS', default='',
                        help='Pefix of the bootstraping results. Default is an empty string, i.e. analysis won\'t be done.')
    result.add_argument('-bsn', '--bs-num', type=int, metavar='NUM', default=50,
                        help='Number of bootstraping. Default is 50')
    result.add_argument('-bsfr', '--bs-fraction', type=float, metavar='PERCENT', default=0.8,
                        help='Percentage of data extracted in each bootstraping. Default is 0.8')

    args = parser.parse_args()
    if len(sys.argv) == 1:  # if cmd input is empty, print help and exit as FEbuilder -h does.
        parser.print_help(sys.stderr)
        sys.exit(1)
    return args


def define_logger(logfile):
    # copied from FEbuilder
    # Create a logger object
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    # Create a handler for writing log messages to a file
    logfile = os.path.join(os.path.abspath('.'), logfile)
    file_handler = logging.FileHandler(logfile)
    file_handler.setLevel(logging.INFO)
    # Create a handler for writing log messages to the console
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    # Add the handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    return logger


def get_final_dg_fepout(fepout):
    """
    Get final dG from the last line of fepout file.
    :param fepout: filename
    :return: dG value
    """
    if type(fepout) == str:
        # return float(open(fepout, 'r').readlines()[-1].strip().split(' ')[18])
        # now 25% faster?
        dq = collections.deque(open(fepout, 'r', encoding='utf-8'))
        last_row = dq.pop()
        try:
            return float(last_row.strip().split(' ')[18])
        except IndexError:  # this line empty?
            last_row = dq.pop()
            return float(last_row.strip().split(' ')[18])
    elif type(fepout) == list:
        return np.array([get_final_dg_fepout(f) for f in fepout])
    else:
        exit('Invalid fepout type!')


## Preprocessing, and tools
def namd_combine_fwbk(u_nk_fwd_lst, u_nk_rev_lst, edit_lambda=True, method=None):
    """
    Sort by lambda value and timesteps; unify column names, especially transform backward lambdas to 1 - lambda
    Avoid numerical precision errors for very small numbers. Nobody will set lambda to 1e-8....

    Combine values in forward and backward DataFrames. For column lambda_i, collect values for lambda_i-1 and lambda_i+1 (except 0 and 1).

    :param u_nk_fwd_lst: list of forward  DataFrames
    :param u_nk_rev_lst: list of backward DataFrames
    :param edit_lambda: edit lambda values for backward data
    :param method: BAR, MBAR or TI
    :return: the combined DataFrame or the result.
    """
    u_nk_fwd = alchemlyb.concat(u_nk_fwd_lst[0])
    u_nk_rev = alchemlyb.concat(u_nk_rev_lst[0])
    u_nk_fwd = u_nk_fwd.sort_index(level=u_nk_fwd.index.names[1:])
    u_nk_fwd.columns = [round(c, 8) for c in u_nk_fwd.columns]
    # reverse lambda values
    if edit_lambda:
        u_nk_rev.columns = [round(1-c, 8) for c in u_nk_rev.columns]  # Update column names
        u_nk_rev.index = u_nk_rev.index.set_levels([round(1-c, 9) for c in u_nk_rev.index.levels[1]], level=1)  # Update fep-lambda indices
        u_nk_rev = u_nk_rev.sort_index(level=u_nk_fwd.index.names[1:]).sort_index(axis=1)  # sort
    # combine the two
    u_nk_fwd.replace(0, np.nan, inplace=True)
    u_nk_fwd[u_nk_fwd.isnull()] = u_nk_rev
    u_nk_fwd.replace(np.nan, 0, inplace=True)
    u_nk = u_nk_fwd.sort_index(level=u_nk_fwd.index.names[1:])
    # sort final dataframe by `fep-lambda` (as opposed to `timestep`)
    if method is None:  # just combine forward/backward
        return u_nk
    else:  # fit. est: estimator
        est = method().fit(u_nk)
        delta_f = to_kcalmol(est.delta_f_)
        return delta_f


def process_sub_diagonal(df):
    """
    :return: the stderr of each window, and the square root of sum of square of these elements.
    """
    # Ensure the DataFrame is square
    if df.shape[0] != df.shape[1]:
        raise ValueError("DataFrame must be square")
    # Calculate the sum of squares of the sub-diagonal elements
    sub_diagonal_elements = np.array([df.iloc[i+1, i] for i in range(df.shape[0]-1)])
    sde = np.insert(sub_diagonal_elements, 0, 0)  # insert 0 for 1st lambda
    norms = np.array([np.linalg.norm(sde[:i]) for i in range(df.shape[0])])
    return sde, norms


def plot_common(xlabel, ylabel, ax=None):
    if ax is None:  # pragma: no cover
        fig, ax = plt.subplots(figsize=(8, 6))
    for dire in ["top", "right", "bottom", "left"]:
        plt.setp(ax.spines[dire], color="#D2B9D3", lw=3, zorder=-2)
    # for dire in ["top", "right"]:
    #     ax.spines[dire].set_color("none")
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.tick_params(axis="both", color="#D2B9D3", labelsize=14, labelcolor='black')
    ax.set_xlabel(xlabel=xlabel, fontsize=16)  # , color="#151B54"
    ax.set_ylabel(ylabel=ylabel, fontsize=16)
    return ax


def plot_dg_fepout(dg_fepout_df):
    ax = plot_common(xlabel="State", ylabel="\u0394G (kcal/mol)")
    sns.set_theme(style="ticks", palette="pastel")
    sns.boxplot(ax=ax, data=dg_fepout_df, x="State", y="dg", hue="Direction", palette=["m", "g"], showfliers=False)
    sns.stripplot(ax=ax, data=dg_fepout_df, x="State", y="dg", hue="Direction", dodge=True, color="black", size=6, jitter=True, legend=False)
    plt.ylabel("\u0394G (kcal/mol)")
    # ax.legend(prop=FP(family='Arial', size=14))  # loses title 'Direction'
    # sns.despine(offset=10, trim=True)
    # plt.show()


def plot_result(dataframe, mode, units=None, final_error=None, ax=None):
    """Plot the forward and backward convergence.

     The columns are: term1, term1_error, term2, term2_error

     Returns
     -------
     matplotlib.axes.Axes
         An axes with the forward and backward convergence drawn.

     Ref
     ----
     The code is taken and modified from
     `alchemlyb <https://github.com/alchemistry/alchemlyb>`_.

    """
    if mode == 'dg-lambda':
        xlabel = r"FEP $\lambda$"
        xticks = np.linspace(0, 1, 11)
    elif mode == 'time-conv':
        xlabel = r"Simulation time per window (ns)"
        xticks = dataframe.index.values
    else:
        raise ValueError("Wrong mode!")

    if units is None:
        units = dataframe.attrs['energy_unit']
    ylabel = r"$\Delta G$ ({})".format(units)
    dataframe = get_unit_converter(units)(dataframe)
    term1 = dataframe.iloc[:, 0].to_numpy()
    term1_error = dataframe.iloc[:, 1].to_numpy()
    term2 = dataframe.iloc[:, 2].to_numpy()
    term2_error = dataframe.iloc[:, 3].to_numpy()

    ax = plot_common(xlabel, ylabel, ax)
    color1 = "#736AFF"
    color2 = "#C11B17"
    # if final_error is None:
    #     final_error = term2_error[-1]
    #
    # if np.isfinite(term2[-1]) and np.isfinite(final_error):
    #     line0 = ax.fill_between(
    #         [0, 1],
    #         term2[-1]-final_error,
    #         term2[-1]+final_error,
    #         color="#D2B9D3",
    #         zorder=1,
    #     )
    line1 = ax.errorbar(
        dataframe.index.values, term1, yerr=term1_error,
        color=color1,
        lw=3,
        zorder=2,
        marker="o",
        mfc="w",
        mew=2.5,
        mec=color1,
        ms=12,
    )
    line2 = ax.errorbar(
        dataframe.index.values, term2, yerr=term2_error,
        color=color2,
        lw=3,
        zorder=3,
        marker="o",
        mfc="w",
        mew=2.5,
        mec=color2,
        ms=12,
    )
    plt.xticks(xticks, [f"{i:.2f}" for i in xticks], fontsize=14)
    plt.yticks(fontsize=14)
    ax.legend(
        (line1[0], line2[0]),
        (dataframe.columns[0], dataframe.columns[2]),
        loc=9,
        prop=FP(size=16),
        frameon=False,
    )
    return ax


class FepoutAnalysis(object):
    def __init__(self, args):
        """
        args (Namespace): Arguments parsed from command line
        filelist (dict): Dictionary containing nested dictionaries for 'bound' and 'free' states,
                             each with lists of file paths for 'forward' and 'backward'.
        Preprocessing: retrieving necessary arguments. Initializes the result folder and data dictionary.
        """
        self.args = args
        self.args.method = methods[args.method]
        if os.path.exists(self.args.result):
            shutil.rmtree(self.args.result)  # TODO: cannot set result folder as . !!
        os.mkdir(self.args.result)
        self.logger = define_logger(os.path.join(self.args.result, self.args.logfile))
        self.filelist, self.num_replicas = self.read_file_list()
        self.lambda_list, self.time_per_window = self.overview_fepout()
        self.result = {
            "deltaf": {},
            "stderr": {},
            "profile": {},
            "ddeltaf": {},
            "time-conv": {},
            "bs-deltaf": {}
        }
        for key in self.result.keys():
            self.result[key]['bound'] = {}
            self.result[key]['free'] = {}

    def run(self):
        self.check_dg()
        self.get_results()
        self.output_results()

    def read_file_list(self):
        """
        Reads and formats the file with paths to the fepout files for both bound and free states.
        Processes only lines with 'forward' or 'backward' keys and 'bound' or 'free' states.
        Skips lines that do not fit these criteria or start with '#'.

        :return: A dictionary containing 'bound' and 'free' states,
              each with 'forward' and 'backward' keys mapped to lists of file paths.
        """
        filelist = {state: {key: [] for key in expected_keys} for state in expected_states}
        with open(self.args.fepout, 'r') as file:
            for line in file:
                # Skip comments and empty lines
                if line.startswith('#') or not line.strip():
                    continue
                state, key, path = line.strip().split(',')
                if state in expected_states and key in expected_keys:
                    filelist[state][key].append(path)
                else:
                    # Optionally, you can print a warning or log skipped lines
                    self.logger.info(f"Skipping line: {line.strip()}")
        # print a summary chart of number of files. Adjust as needed for longer names
        state_width = 7
        forward_width = 15
        backward_width = 15
        self.logger.info("Number of replicas should be the same:")
        self.logger.info(f"| {'State'.ljust(state_width)} | {'Forward Files'.ljust(forward_width)} | {'Backward Files'.ljust(backward_width)} |")
        self.logger.info(f"| {'-'*state_width} | {'-'*forward_width} | {'-'*backward_width} |")
        for state in expected_states:
            forward_count = len(filelist[state]['forward'])
            backward_count = len(filelist[state]['backward'])
            self.logger.info(f"| {state.capitalize().ljust(state_width)} | {str(forward_count).ljust(forward_width)} | {str(backward_count).ljust(backward_width)} |")
        self.logger.info('')
        return filelist, len(filelist[state]['forward'])

    def overview_fepout(self):
        """
        Take one file and retrieve the lambdas.
        :return: the list of lambdas.
        """
        data = apn.extract_u_nk(self.filelist['bound']['forward'][0], self.args.temperature)
        lambda_list = data.columns.tolist()
        time_idx = data.index.get_level_values(0)
        time_per_window = (max(time_idx) - min(time_idx))/1000000*self.args.dt  # unit: nanosecond
        return lambda_list, time_per_window

    def check_dg(self):
        paths = []
        for state, directions in self.filelist.items():
            for direction, files in directions.items():
                for file in files:
                    paths.append({'State': state, 'Direction': direction, 'dg': get_final_dg_fepout(file)})
        dg_fepout_df = pd.DataFrame(paths)
        dg_fepout_df['dg'] = dg_fepout_df.apply(lambda row: row['dg']*(-1) if row['Direction'] == 'backward' else row['dg'], axis=1)
        plot_dg_fepout(dg_fepout_df)
        plt.savefig(os.path.join(self.args.result, self.args.fodgfile+'.png'))
        dg_fepout_df.to_csv(os.path.join(str(self.args.result), self.args.fodgfile+'.csv'))
        self.logger.info("\u0394G from .fepout files plotted. Check outliers and remove them from input file list!\n")

    # fitting for free energy change
    def get_results(self):
        """
        Get the results from the fepout files, handling multiple file paths for each key/state.

        :return: No, just update the result dictionary.
        """
        # Iterate over states ('bound' and 'free')
        for state in self.filelist.keys():
            fpaths = self.filelist[state]['forward']
            bpaths = self.filelist[state]['backward']
            self.update_results(state, fpaths, bpaths)
            if not self.args.no_decompose:
                self.update_decomposition(state, fpaths, bpaths, 'elec')
                self.update_decomposition(state, fpaths, bpaths, 'vdw')
                self.logger.info(f'Done {state} decomposition')

    def update_results(self, state, fpaths, bpaths):
        """
        Update the results for each state.

        Parameters:
        :param state (str): bound or free.
        :param fpath (list): List of file paths for forward calculations.
        :param bpath (list): List of file paths for backward calculations.
        """
        # Total
        f_df = [apn.extract_u_nk(fp, self.args.temperature, self.lambda_list) for fp in fpaths]
        b_df = [apn.extract_u_nk(bp, self.args.temperature, self.lambda_list) for bp in bpaths]
        u_nk = namd_combine_fwbk(f_df, b_df, method=None)
        if self.args.convfile == '':
            est = self.fit_bidirection(u_nk)
            self.update_dict(est, state, 'dG')
        else:
            self.logger.info(f'Doing time convergence analysis for {state}...')
            est_list = self.fit_bidirection(u_nk, time_conv=True)
            self.update_dict(est_list[-1], state, 'dG')
            self.update_time_convergence(est_list, state)
        self.logger.info(f'Done {state} dG')

        # boostraping
        if self.args.bsfile != '':
            self.result["bs-deltaf"][state] = self.fit_bootstrap(u_nk)
            self.logger.info(f'Done {state} bootstrapping')

    def update_decomposition(self, state, fpaths, bpaths, component):
        # separate this to avoid memory overflow...
        # Electrostatics and van der Waals components
        df_fwd = [apn.extract_u_nk(fp, self.args.temperature, self.lambda_list, component) for fp in fpaths]
        df_rev = [apn.extract_u_nk(bp, self.args.temperature, self.lambda_list, component) for bp in bpaths]
        u_nk_component = namd_combine_fwbk(df_fwd, df_rev, method=None)
        est = self.fit_bidirection(u_nk_component)
        self.update_dict(est, state, component)

    def fit_bidirection(self, u_nk, time_conv=False):
        """
        Do the fitting.

        :param u_nk: combined forward/backward data returned by `namd_combine_fwbk`
        :param component: dG, elec or vdw.
        :return: an estimator object, or a list of them for time convergence analysis
        """
        if not time_conv:
            est = self.args.method().fit(u_nk)
            return est
        else:
            est_list = []
            num_sampling_per_window = len(u_nk[0]) / len(self.lambda_list) / self.num_replicas
            for i in range(1, self.args.conv_num+1):
                subsampling = num_sampling_per_window/self.args.conv_num*i
                u_nk_cut_list = [u_nk.iloc[int(j * num_sampling_per_window * self.num_replicas):
                                           int(j * num_sampling_per_window * self.num_replicas + subsampling * self.num_replicas)]
                                 for j in range(len(self.lambda_list))]
                u_nk_cut = alchemlyb.concat(u_nk_cut_list)
                est_list.append(self.args.method().fit(u_nk_cut))
                # self.logger.info(f'Done {i/self.args.conv_num*100:.1f}%')
            return est_list

    def fit_bootstrap(self, u_nk):
        results = []
        for i in range(self.args.bs_num):
            # subsample a fraction of the data randomly, in each value of fep-lambda
            # remove the index added, and sort by time
            u_nk_subsample = (u_nk.groupby('fep-lambda').apply(
                lambda group: group.sample(frac=self.args.bs_fraction)).
                reset_index(level=0, drop=True).sort_index(level=u_nk.index.names[1:]))
            est = self.args.method().fit(u_nk_subsample)
            dg = to_kcalmol(est.delta_f_)[1.0][0.0]
            results.append(dg)
        return np.array(results, dtype=float)

    def update_dict(self, est, state, component):
        """
        Actually update the results.

        :param est: returned by `fit_bidirection`.
        :param state:
        :param component:
        """
        delta_f = to_kcalmol(est.delta_f_)
        ddeltaf = to_kcalmol(est.d_delta_f_)
        sde, norms = process_sub_diagonal(ddeltaf)
        self.result["deltaf"][state][component] = delta_f[1.0][0.0]
        self.result["stderr"][state][component] = norms[-1]
        self.result["profile"][state][component] = delta_f[0.0].values
        self.result["ddeltaf"][state][component] = norms

    def update_time_convergence(self, est_list, state):
        deltaf_list = []
        stderr_list = []
        for est in est_list:
            deltaf_list.append(to_kcalmol(est.delta_f_)[1.0][0.0])
            ddeltaf = to_kcalmol(est.d_delta_f_)
            sde, norms = process_sub_diagonal(ddeltaf)
            stderr_list.append(norms[-1])
        # make the data a DataFrame
        self.result['time-conv'][state] = pd.DataFrame({
            "deltaf": deltaf_list,
            "stderr": stderr_list
        }, index=np.linspace(1/self.args.conv_num, 1, self.args.conv_num)*self.num_replicas*self.time_per_window)

    # output
    def output_results(self):
        """
        Print the results
        """
        # summary of energy
        self.logger.info('\nResult:')
        self.logger.info("Forward \u0394\u0394G is {:9.3f} ± {:.3f} kcal/mol".format(
            self.result["deltaf"]['bound']['dG']-self.result["deltaf"]['free']['dG'],
            np.linalg.norm([self.result["stderr"]['bound']['dG'], self.result["stderr"]['free']['dG']])))
        if self.args.no_decompose:
            self.logger.info(f'Disabled decomposition')
        else:
            for component in expected_components:
                self.logger.info("\u0394G_{:5s}    is {:9.3f} ± {:.3f} kcal/mol".format(component,
                            self.result["deltaf"]['bound'][component]-self.result["deltaf"]['free'][component],
                            np.linalg.norm([self.result["stderr"]['bound'][component], self.result["stderr"]['free'][component]])))


        # dG-lambda plot
        dg_lambda_df = pd.DataFrame({
            "Bound": self.result["profile"]['bound']['dG'],
            "Bound_error": self.result["ddeltaf"]['bound']['dG'],
            "Free": self.result["profile"]['free']['dG'],
            "Free_error": self.result["ddeltaf"]['free']['dG']
        }, index=self.lambda_list)
        dg_lambda_df.attrs = {'temperature': self.args.temperature, 'energy_unit': 'kcal/mol'}
        ax = plot_result(dg_lambda_df, mode='dg-lambda')
        ax.figure.savefig(os.path.join(self.args.result, self.args.dglfile+'.png'))
        dg_lambda_df.to_csv(os.path.join(str(self.args.result), self.args.dglfile+'.csv'))

        # time convergence (bound/free)
        if not self.args.convfile == '':
            time_conv_df = pd.concat([self.result['time-conv']['bound'], self.result['time-conv']['free']], axis=1)
            time_conv_df.columns = ['Bound', 'Bound_error', 'Free', 'Free_error']
            time_conv_df.attrs = {'temperature': self.args.temperature, 'energy_unit': 'kcal/mol'}
            ax = plot_result(time_conv_df, mode='time-conv')
            ax.figure.savefig(os.path.join(self.args.result, self.args.convfile+'.png'))
            time_conv_df.to_csv(os.path.join(str(self.args.result), self.args.convfile+'.csv'))

        # bootstraping (bound/free)
        if self.args.bsfile != '':
            self.logger.info("Bootstrapping results:")
            string = ''
            for state in expected_states:
                self.logger.info("\u0394G_{:5s}    is {:9.3f} ± {:.3f} kcal/mol".format(state,
                            self.result["deltaf"][state]['dG'], np.std(self.result["bs-deltaf"][state])))
                str_list = list(map(str, self.result["bs-deltaf"][state]))
                string += "\u0394G_{:5s},{:s}\n".format(state, ','.join(str_list))
            _ = open(os.path.join(self.args.result, self.args.bsfile+'.csv'), 'w').write(string)


def main():
    args = parse_arguments(__version__)
    fa = FepoutAnalysis(args)
    fa.run()



if __name__ == "__main__":
    main()


# %% debug
# if __name__ == "__main__":
#     # path = '/data/work/make_hybrid_top/tests/Lei2Guanqiao'
#     # path = '/data/work/make_hybrid_top/openforcefield-protein-ligand-benchmark/hif2a/FEP/prod1'
#     # path = 'E:\\GitHub-repo\\make_hybrid_top\\test\\Lei2Guanqiao'
#     # path = 'E:\\GitHub-repo\\make_hybrid_top\\test\\openforcefield-protein-ligand-benchmark\\hif2a\\FEP\\prod1'
#     # path = 'E:\\GitHub-repo\\make_hybrid_top\\test\\fm2gxf\\B44\\SEIWDRIDF2AEIAARAAA'
#     path = 'E:\\GitHub-repo\\make_hybrid_top\\test\\fm2gxf\\BA\\GILGFVFTL2KLWAQCVQL'
#     os.chdir(path)
#     data = {
#         # 'fepout': os.path.join(path, 'D3N-N3D.csv'),
#         # 'fepout': os.path.join(path, '42-10-10-42.csv'),
#         'fepout': os.path.join(path, 'forward-backward.csv'),
#         'method': 'BAR',
#         'temperature': 310,
#         # 'result': 'D3N-N3D',
#         # 'result': '42-10-10-42',
#         'result': 'result',
#         'logfile': 'analysis.log',
#         "fodgfile": 'dg_fepout',
#         "no_decompose": False,  # False True
#         "dglfile": 'dg-lambda',
#         "convfile": '',  # time-conv
#         "conv_num": 10,
#         "dt": 1.0,
#         "bsfile": "",  # bootstrap
#         "bs_num": 50,
#         "bs_fraction": 0.8
#     }
#     args = argparse.Namespace(**data)
#     fa = FepoutAnalysis(args)
#     fa.run()