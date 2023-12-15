#!/usr/bin/env python

######################################################################
# Alchemical Analysis with ChatGPT: A tool for analyzing alchemical free energy calculations using ChatGPT
# Copyright 2023 Zhejiang University and the Authors
#
# Authors: KCC
#
# This code is a modified version of the Alchemical Analysis tool developed by Pavel Klimovich, Michael Shirts, and David Mobley.
# It has been adapted to incorporate ChatGPT, a language model developed by OpenAI, for enhanced analysis capabilities.
#
# This code is distributed under the GNU Lesser General Public License, version 2.1 or any later version.
#
# This code is provided "as is" without warranty of any kind, expressed or implied, including but not limited to the
# warranties of merchantability, fitness for a particular purpose, and non-infringement. In no event shall the authors or
# copyright holders be liable for any claim, damages, or other liability, whether in an action of contract, tort, or
# otherwise, arising from, out of, or in connection with this code or the use or other dealings in this code.
#
# If you use this code, please acknowledge the original authors by citing the relevant publications.
# You should have received a copy of the GNU Lesser General Public License along with this code.
# If not, see <http://www.gnu.org/licenses/>.
######################################################################

#===================================================================================================
# IMPORTS
#===================================================================================================

## Not a built-in module. Will be called from main, whenever needed. ##
## import pymbar      Multistate Bennett Acceptance Ratio estimator. ##

import sys
import numpy
import pickle                       # for full-precision data storage
from   optparse import OptionParser # for parsing command-line options
import os                           # for os interface
import time as ttt_time             # for timing
import pdb                          # for debugging
from utils.zeroxvg import *

#===================================================================================================
# INPUT OPTIONS
#===================================================================================================

parser = OptionParser()
parser.add_option('-b', '--fraction', dest = 'fraction', help = 'The fraction of the energy file will be used to calculate the statistics. Default: 1.0',  default = 1.0, type = float)
parser.add_option('-d', '--dir', dest = 'datafile_directory', help = 'Directory in which data files are stored. Default: Current directory.', default = '.')
parser.add_option('-e', '--backward', dest = 'backward', help = 'Extract the energy data from the backward direction. Default: False', default = False, action = 'store_true')
parser.add_option('-f', '--forwrev', dest = 'bForwrev', help = 'Plot the free energy change as a function of time in both directions, with the specified number of points in the time plot. The number of time points (an integer) must be provided. Default: 0', default = 0, type=int)
parser.add_option('-g', '--breakdown', dest = 'breakdown', help = 'Plot the free energy differences evaluated for each pair of adjacent states for all methods, including the dH/dlambda curve for TI. Default: False.', default = False, action = 'store_true')
parser.add_option('-i', '--threshold', dest = 'uncorr_threshold', help = 'Proceed with correlated samples if the number of uncorrelated samples is found to be less than this number. If 0 is given, the time series analysis will not be performed at all. Default: 50.', default = 50, type=int)
parser.add_option('-j', '--resultfilename', dest = 'resultfilename', help = 'custom defined result filename prefix. Default: results', default = 'results')
parser.add_option('-k', '--koff', dest = 'bSkipLambdaIndex', help = 'Give a string of lambda indices separated by \'-\' and they will be removed from the analysis. (Another approach is to have only the files of interest present in the directory). Default: None.', default = '')
parser.add_option('-m', '--methods', dest = 'methods', help = 'A list of the methods to esitimate the free energy with. Default: [TI, TI-CUBIC, DEXP, IEXP, BAR, MBAR]. To add/remove methods to the above list provide a string formed of the method strings preceded with +/-. For example, \'-ti_cubic+gdel\' will turn methods into [TI, DEXP, IEXP, BAR, MBAR, GDEL]. \'ti_cubic+gdel\', on the other hand, will call [TI-CUBIC, GDEL]. \'all\' calls the full list of supported methods [TI, TI-CUBIC, DEXP, IEXP, GINS, GDEL, BAR, UBAR, RBAR, MBAR].', default = '')
parser.add_option('-n', '--uncorr', dest = 'uncorr', help = 'The observable to be used for the autocorrelation analysis; either \'dhdl_all\' (obtained as a sum over all energy components) or \'dhdl\' (obtained as a sum over those energy components that are changing; default) or \'dE\'. In the latter case the energy differences dE_{i,i+1} (dE_{i,i-1} for the last lambda) are used.', default = 'dhdl')
parser.add_option('-o', '--out', dest = 'output_directory', help = 'Directory in which the output files produced by this script will be stored. Default: Same as datafile_directory.', default = '')
parser.add_option('-p', '--prefix', dest = 'prefix', help = 'Prefix for datafile sets, i.e.\'dhdl\' (default).', default = 'dhdl')
parser.add_option('-q', '--suffix', dest = 'suffix', help = 'Suffix for datafile sets, i.e. \'xvg\' (default).', default = 'xvg')
parser.add_option('-r', '--decimal', dest = 'decimal', help = 'The number of decimal places the free energies are to be reported with. No worries, this is for the text output only; the full-precision data will be stored in \'results.pickle\'. Default: 3.', default = 3, type=int)
parser.add_option('-s', '--skiptime', dest = 'equiltime', help = 'Discard data prior to this specified time as \'equilibration\' data. Units picoseconds. Default: 0 ps.', default = 0, type=float)
parser.add_option('-t', '--temperature', dest = 'temperature', help = "Temperature in K. Default: 298 K.", default = 298, type=float)
parser.add_option('-u', '--units', dest = 'units', help = 'Units to report energies: \'kJ\', \'kcal\', and \'kBT\'. Default: \'kJ\'', default = 'kJ')
parser.add_option('-v', '--verbose', dest = 'verbose', help = 'Verbose option. Default: False.', default = False, action = 'store_true')
parser.add_option('-w', '--overlap', dest = 'overlap', help = 'Print out and plot the overlap matrix. Default: False.', default = False, action = 'store_true')
parser.add_option('-x', '--ignoreWL', dest = 'bIgnoreWL', help = 'Do not check whether the WL weights are equilibrated. No log file needed as an accompanying input.', default = False, action = 'store_true')
parser.add_option('-y', '--tolerance', dest = 'relative_tolerance', help = "Convergence criterion for the energy estimates with BAR and MBAR. Default: 1e-10.", default = 1e-10, type=float)
parser.add_option('-z', '--initialize', dest = 'init_with', help = 'The initial MBAR free energy guess; either \'BAR\' or \'zeros\'. Default: \'BAR\'.', default = 'BAR')

import pymbar

#===================================================================================================
# FUNCTIONS: Miscellanea.
#===================================================================================================
def check_units_and_more(units):
    kB = 1.3806488*6.02214129/1000.0  # Boltzmann's constant (kJ/mol/K).
    beta = 1./(kB*P.temperature)

    if units == 'kJ':
        beta_report = beta
        units = '(kJ/mol)'
    elif units == 'kcal':
        beta_report = beta * 4.184
        units = '(kcal/mol)'
    elif units == 'kBT':
        beta_report = 1
        units = '(k_BT)'
    else:
        raise ValueError(f"I don't understand the unit type '{units}': the only options 'kJ', 'kcal', and 'kBT'")

    P.output_directory = P.output_directory or P.datafile_directory

    if P.overlap and 'MBAR' not in P.methods:
        raise ValueError("MBAR is not in 'methods'; can't plot the overlap matrix.")
    
    return units, beta, beta_report

from datetime import timedelta

def timeStatistics(stime):
    elapsed_time = timedelta(seconds=int(ttt_time.time() - stime))
    th, remainder = divmod(elapsed_time.seconds, 3600)
    tm, ts = divmod(remainder, 60)
    current_time = ttt_time.asctime()
    return th, tm, f'{ts:.2f}', current_time

#===================================================================================================
# FUNCTIONS: The autocorrelation analysis.
#===================================================================================================

def uncorrelate(sta, fin, do_dhdl=False):
    if not P.uncorr_threshold:
        return (dhdlt, nsnapshots, None) if P.software.title()=='Sire' else (dhdlt, nsnapshots, u_klt)

    u_kln = numpy.zeros([K,K,max(fin-sta)], numpy.float64)
    N_k = numpy.zeros(K, int)
    g = numpy.zeros(K,float)
    if do_dhdl:
        dhdl = numpy.zeros([K,n_components,max(fin-sta)], float)
        print(f"\n\nNumber of correlated and uncorrelated samples:\n\n{'State':6s} {'N':12s} {'N_k':12s} {'N/N_k':12s}\n")

    for k in range(K):
        lastl = k
        for l in range(K):
            if numpy.array_equal(lv[k],lv[l]):
                lastl = l
        dhdl_sum = numpy.sum(dhdlt[k, lchange[lastl], sta[k]:fin[k]], axis=0)

        g[k] = 1 if not numpy.any(dhdl_sum) else pymbar.timeseries.statisticalInefficiency(dhdl_sum)

        indices = sta[k] + numpy.array(pymbar.timeseries.subsampleCorrelatedData(dhdl_sum, g=g[k]))
        N_uncorr = len(indices)

        if N_uncorr < P.uncorr_threshold:
            if do_dhdl:
                print(f"WARNING: Only {N_uncorr} uncorrelated samples found at lambda number {k}; proceeding with analysis using correlated samples...")
            indices = sta[k] + numpy.arange(len(dhdl_sum))
            N = len(indices)
        else:
            N = N_uncorr
        N_k[k] = N
        if u_klt is not None:
            for l in range(K):
                u_kln[k,l,0:N] = u_klt[k,l,indices]
        if do_dhdl:
            print(f"{k:6s} {N_uncorr:12s} {N_k[k]:12s} {g[k]:12.2f}")
            for n in range(n_components):
                dhdl[k,n,0:N] = dhdlt[k,n,indices]

    return (dhdl, N_k, u_kln) if do_dhdl else (N_k, u_kln)

#===================================================================================================
# FUNCTIONS: The MBAR workhorse.
#===================================================================================================

def estimate_with_MBAR(u_kln, N_k, reltol, regular_estimate=False):
    """Estimate the free energy using the Multistate Bennett Acceptance Ratio (MBAR) method."""

    def plotOverlapMatrix(O):
		"""Plots the probability of observing a sample from state i (row) in state j (column).
		For convenience, the neigboring state cells are fringed in bold."""
		max_prob = O.max()
		fig = pl.figure(figsize=(K/2.,K/2.))
		fig.add_subplot(111, frameon=False, xticks=[], yticks=[])

		for i in range(K):
			if i!=0:
				pl.axvline(x=i, ls='-', lw=0.5, color='k', alpha=0.25)
				pl.axhline(y=i, ls='-', lw=0.5, color='k', alpha=0.25)
			for j in range(K):
				if O[j,i] < 0.005:
					ii = ''
				elif O[j,i] > 0.995:
				    ii = '1.00'
				else:
				    ii = ("%.2f" % O[j,i])[1:]
				alf = O[j,i]/max_prob
				pl.fill_between([i,i+1], [K-j,K-j], [K-(j+1),K-(j+1)], color='k', alpha=alf)
				pl.annotate(ii, xy=(i,j), xytext=(i+0.5,K-(j+0.5)), size=8, textcoords='data', va='center', ha='center', color=('k' if alf < 0.5 else 'w'))

		if P.bSkipLambdaIndex:
			ks = [int(l) for l in P.bSkipLambdaIndex.split('-')]
			ks = numpy.delete(numpy.arange(K+len(ks)), ks)
		else:
			ks = range(K)
		for i in range(K):
			pl.annotate(ks[i], xy=(i+0.5, 1), xytext=(i+0.5, K+0.5), size=10, textcoords=('data', 'data'), va='center', ha='center', color='k')
			pl.annotate(ks[i], xy=(-0.5, K-(j+0.5)), xytext=(-0.5, K-(i+0.5)), size=10, textcoords=('data', 'data'), va='center', ha='center', color='k')
		pl.annotate('$\lambda$', xy=(-0.5, K-(j+0.5)), xytext=(-0.5, K+0.5), size=10, textcoords=('data', 'data'), va='center', ha='center', color='k')
		pl.plot([0,K], [0,0], 'k-', lw=4.0, solid_capstyle='butt')
		pl.plot([K,K], [0,K], 'k-', lw=4.0, solid_capstyle='butt')
		pl.plot([0,0], [0,K], 'k-', lw=2.0, solid_capstyle='butt')
		pl.plot([0,K], [K,K], 'k-', lw=2.0, solid_capstyle='butt')

		cx = sorted(2*range(K+1))
		cy = sorted(2*range(K+1), reverse=True)
		pl.plot(cx[2:-1], cy[1:-2], 'k-', lw=2.0)
		pl.plot(numpy.array(cx[2:-3])+1, cy[1:-4], 'k-', lw=2.0)
		pl.plot(cx[1:-2], numpy.array(cy[:-3])-1, 'k-', lw=2.0)
		pl.plot(cx[1:-4], numpy.array(cy[:-5])-2, 'k-', lw=2.0)

		pl.xlim(-1, K)
		pl.ylim(0, K+1)
		pl.savefig(os.path.join(P.output_directory, 'O_MBAR.pdf'), bbox_inches='tight', pad_inches=0.0)
		pl.close(fig)
		return

    if regular_estimate:
        print("\nEstimating the free energy change with MBAR...")
    MBAR = pymbar.mbar.MBAR(u_kln, N_k, verbose=P.verbose, relative_tolerance=reltol, initialize=P.init_with)

    Deltaf_ij, dDeltaf_ij, theta_ij = MBAR.getFreeEnergyDifferences(uncertainty_method='svd-ew', return_theta=True)

    if P.verbose:
        print(f"Matrix of free energy differences\nDeltaf_ij:\n{Deltaf_ij}\ndDeltaf_ij:\n{dDeltaf_ij}")

    if regular_estimate and P.overlap:
        print("The overlap matrix is...")
        O = MBAR.computeOverlap()[2]
        for k in range(K):
            line = ' '.join(f'{O[k, l]:5.2f}' for l in range(K))
            print(line)
        plotOverlapMatrix(O)
        print("\nFor a nicer figure look at 'O_MBAR.pdf'")
        return Deltaf_ij, dDeltaf_ij

    return Deltaf_ij[0,K-1]/P.beta_report, dDeltaf_ij[0,K-1]/P.beta_report

#===================================================================================================
# FUNCTIONS: This one estimates dF and ddF for all pairs of adjacent states and stores them.
#===================================================================================================

def estimate_pairs():
    """Estimate the free energy differences and their uncertainties using various methods."""
    print(f"Estimating the free energy change with {', '.join(P.methods)}")

    df_allk = list()
    ddf_allk = list()

    for k in range(K-1):
        df = dict()
        ddf = dict()

        for name in P.methods:
            if name == 'TI':
	            #===================================================================================================
	            # Estimate free energy difference with TI; interpolating with the trapezoidal rule.
	            #===================================================================================================
                df['TI'] = 0.5 * numpy.dot(dlam[k], (ave_dhdl[k] + ave_dhdl[k+1]))
                ddf['TI'] = 0.5 * numpy.sqrt(numpy.dot(dlam[k]**2, std_dhdl[k]**2 + std_dhdl[k+1]**2))

            elif name == 'TI-CUBIC':
            	#===================================================================================================
	            # Estimate free energy difference with TI; interpolating with the natural cubic splines.
	            #===================================================================================================
                df['TI-CUBIC'], ddf['TI-CUBIC'] = 0, 0
	            for j in range(n_components):
	               if dlam[k,j] > 0:
	                  lj = lchange[:,j]
	                  df['TI-CUBIC'] += numpy.dot(cubspl[j].wk[mapl[k,j]],ave_dhdl[lj,j])
	                  ddf['TI-CUBIC'] += numpy.dot(cubspl[j].wk[mapl[k,j]]**2,std_dhdl[lj,j]**2)
	            ddf['TI-CUBIC'] = numpy.sqrt(ddf['TI-CUBIC'])
                
            elif name in ['DEXP', 'GDEL']:
                w_F = u_kln[k, k+1, 0:N_k[k]] - u_kln[k, k, 0:N_k[k]]
                if name == 'DEXP':
                	#===================================================================================================
		            # Estimate free energy difference with Forward-direction EXP (in this case, Deletion from solvent).
		            #===================================================================================================
                    df['DEXP'], ddf['DEXP'] = pymbar.exp.EXP(w_F)
                elif name == 'GDEL':
                	#===================================================================================================
		            # Estimate free energy difference with a Gaussian estimate of EXP (in this case, deletion from solvent)
		            #===================================================================================================
                    df['GDEL'], ddf['GDEL'] = pymbar.exp.EXPGauss(w_F)

            elif name in ['IEXP', 'GINS', 'BAR', 'UBAR', 'RBAR']:
            	w_F = u_kln[k, k+1, 0:N_k[k]] - u_kln[k, k, 0:N_k[k]]
                w_R = u_kln[k+1, k, 0:N_k[k+1]] - u_kln[k+1, k+1, 0:N_k[k+1]]
                if name == 'IEXP':
                	#===================================================================================================
		            # Estimate free energy difference with Reverse-direction EXP (in this case, insertion into solvent).
		            #===================================================================================================
                    rdf, rddf = pymbar.exp.EXP(w_R)
                    df['IEXP'], ddf['IEXP'] = (-rdf, rddf)
                elif name == 'GINS':
                	#===================================================================================================
		            # Estimate free energy difference with a Gaussian estimate of EXP (in this case, insertion into solvent)
		            #===================================================================================================
                    rdf, rddf = pymbar.exp.EXPGauss(w_R)
                    df['GINS'], ddf['GINS'] = (-rdf, rddf)
			    elif name == 'BAR':
			    	#===================================================================================================
		            # Estimate free energy difference with BAR; use w_F and w_R computed above.
		            #===================================================================================================
			        df['BAR'], ddf['BAR'] = pymbar.bar.BAR(w_F, w_R, relative_tolerance=P.relative_tolerance, verbose = P.verbose)
			    elif name == 'UBAR':
			    	#===================================================================================================
		            # Estimate free energy difference with unoptimized BAR -- assume dF is zero, and just do one evaluation
		            #===================================================================================================
			        df['UBAR'], ddf['UBAR'] = pymbar.bar.BAR(w_F, w_R, verbose = P.verbose,iterated_solution=False)
			    elif name == 'RBAR':
			    	#===================================================================================================
		            # Estimate free energy difference with Unoptimized BAR over range of free energy values, and choose the one
		            # that is self consistently best.
		            #===================================================================================================
			        min_diff = 1E6
			        best_udf = 0
			        for trial_udf in range(-10,10,1):
			            udf, uddf = pymbar.bar.BAR(w_F, w_R, DeltaF=trial_udf, iterated_solution=False, verbose=P.verbose)
			            diff = numpy.abs(udf - trial_udf)
			            if (diff < min_diff):
			                best_udf = udf
			                best_uddf = uddf
			                min_diff = diff
			        df['RBAR'], ddf['RBAR'] = (best_udf, best_uddf)
			                
		    if name == 'MBAR':
	            #===================================================================================================
	            # Store the MBAR free energy difference (already estimated above) properly, i.e. by state.
	            #===================================================================================================     
		        df['MBAR'], ddf['MBAR'] =  Deltaf_ij[k, k+1], numpy.nan_to_num(dDeltaf_ij[k, k+1])

        df_allk.append(df)
        ddf_allk.append(ddf)

	return df_allk, ddf_allk

#===================================================================================================
# FUNCTIONS: All done with calculations; summarize and print stats.
#===================================================================================================

def total_energies():

    # Count up the charging states.
    start_coul, end_coul, start_vdw, end_vdw = 0, 0, 0, 0
    for lv_n in ['vdw','coul','fep']:
        if lv_n in P.lv_names:
            _ndx_char = P.lv_names.index(lv_n)
            lv_char = lv[:, _ndx_char]
            if not np.all(lv_char == lv_char[0]):
                if lv_n in ['vdw', 'fep']:
                    end_vdw = end_coul = np.sum(lv_char != 1)
                    ndx_char = _ndx_char
                elif lv_n == 'coul':
                    end_coul = np.sum(lv_char != 1)
                    ndx_char = _ndx_char

    # Determine where Coulomb section is, if it is present.
    if end_coul > end_vdw:
        start_coul = end_vdw
        start_vdw = 0
    elif start_coul == end_coul and P.verbose:
        print("No Coulomb transformation present.")
    else:
        start_coul = 0
        start_vdw = end_coul

    segments = ['Coulomb', 'vdWaals', 'TOTAL']
    segment_starts = [start_coul, start_vdw, 0]
    segment_ends = [min(end_coul, K-1), min(end_vdw, K-1), K-1]
    dFs, ddFs = [], []

    # Perform the energy segmentation; be pedantic about the TI cumulative ddF's (see Section 3.1 of the paper).
    for segment, seg_start, seg_end in zip(segments, segment_starts, segment_ends):
        dF = {name: 0 for name in P.methods}
        ddF = {name: 0 for name in P.methods}

        for name in P.methods:
            if name.startswith('TI'):
                for k in range(seg_start, seg_end):
                    dF[name] += df_allk[k][name]

                jlist = [ndx_char] if segment == 'Coulomb' and end_coul > 0 else []
                if segment == 'vdWaals':
                    jlist = []
                elif segment == 'TOTAL':
                    jlist = list(range(n_components))

                for j in jlist:
                    lj = lchange[:,j]
                    if not np.all(lj == False):  # handle the all-zero lv column
                        if name == 'TI-CUBIC':
                            ddF[name] += np.dot((cubspl[j].wsum)**2, std_dhdl[lj,j]**2)
                        elif name == 'TI':
                            h = np.trim_zeros(dlam[:,j])
                            wsum = 0.5 * (np.append(h, 0) + np.append(0, h))
                            ddF[name] += np.dot(wsum**2, std_dhdl[lj,j]**2)
 
                ddF[name] = np.sqrt(ddF[name])

            elif name == 'MBAR':
                dF[name] = Deltaf_ij[seg_start,seg_end]
                ddF[name] = dDeltaf_ij[seg_start,seg_end]

            else:
                for k in range(seg_start, seg_end):
                    dF[name] += df_allk[k][name]
                    ddF[name] += (ddf_allk[k][name])**2
                ddF[name] = np.sqrt(ddF[name])

        dFs.append(dF)
        ddFs.append(ddF)

    for name in P.methods:  # 'vdWaals' = 'TOTAL' - 'Coulomb'
        ddFs[1][name] = np.sqrt(ddFs[2][name]**2 - ddFs[0][name]**2)

    def print_line(text, data_format, d1=None, d2=None):
        """Fills out the results table linewise."""
        output = [text]
        for name in P.methods:
            if d1 is not None and d2 is not None:
                output.append(data_format % (d1[name] / P.beta_report, d2[name] / P.beta_report))
            elif d1 == 'name':
                output.append(data_format % (name, P.units))
            else:
                output.append(data_format)
        print(' '.join(output))
        out_text.append(' '.join(output) + '\n')

    d = P.decimal
    data_format = ('X%d.%df  +-  X%d.%df' % (d+7, d, d+2, d)).replace('X', '%')
    names_format = ('X%ds X-%ds' % (d+6, d+8)).replace('X', '%')
    out_text = []
    print_line(12*'-', names_format, 'plain')
    print_line('%-12s' % '   States', names_format, 'name')
    print_line(12*'-', names_format, 'plain')
    for k in range(K-1):
        print_line('%4d -- %-4d' % (k, k+1), data_format, df_allk[k], ddf_allk[k])
    print_line(12*'-', names_format, 'plain')
    if len(P.lv_names) > 1:
        for i in range(len(segments)):
            print_line('%9s:  ' % segments[i], data_format, dFs[i], ddFs[i])
        print('{:^40}'.format("\n".join([
            "A remark on the energy components interpretation:",
            "'vdWaals' is computed as 'TOTAL' - 'Coulomb', where",
            "'Coulomb' is found as the free energy change between",
            "the states defined by the lambda vectors (0,0,...,0)",
            "and (1,0,...,0), the only varying vector component",
            "being either 'coul-lambda' or 'fep-lambda'."])))
    else:
        print_line('%9s:  ' % segments[-1], data_format, dFs[-1], ddFs[-1])
    
    # Store results.
    with open(os.path.join(P.output_directory, f'{P.resultfilename}.txt'), 'w') as outfile:
        outfile.write('# Command line was: %s\n\n' % ' '.join(sys.argv))
        outfile.writelines(out_text)

    P.datafile_directory = os.getcwd()
    P.when_analyzed = time.asctime()
    P.ddf_allk = ddf_allk
    P.df_allk  = df_allk
    P.ddFs     = ddFs
    P.dFs      = dFs

    with open(os.path.join(P.output_directory, f'{P.resultfilename}.pickle'), 'wb') as outfile:
        pickle.dump(P, outfile)

    print('\n' + w * '*')
    for i in [" The above table has been stored in ", f" {P.output_directory}/{P.resultfilename}.txt ",
              " while the full-precision data ", " (along with the simulation profile) in ", 
              f" {P.output_directory}/{P.resultfilename}.pickle "]:
        print('{:^40}'.format(i))
    print(w * '*')

    return

#===================================================================================================
# FUNCTIONS: Free energy change breakdown (into lambda-pair dFs). Called by the -g flag.
#===================================================================================================

def plot_dFvsLambda():
	import matplotlib.pyplot as plt
	import numpy as np

	def plot_dF_vs_lambda1():
	    """Plots the free energy differences evaluated for each pair of adjacent states for all methods."""
	    x = np.arange(len(df_allk))
	    fig_size = (8, 6) if x[-1] < 8 else (len(x), 6)
	    fig, ax = plt.subplots(figsize=fig_size)

	    width = 1./(len(P.methods) + 1)
	    elw = 30 * width

	    colors = {'TI':'#C45AEC', 'TI-CUBIC':'#33CC33', 'DEXP':'#F87431', 'IEXP':'#FF3030', 'GINS':'#EAC117', 
	              'GDEL':'#347235', 'BAR':'#6698FF', 'UBAR':'#817339', 'RBAR':'#C11B17', 'MBAR':'#F9B7FF'}
	    lines = []

	    for idx, name in enumerate(P.methods):
	        y = [df_allk[i][name] / P.beta_report for i in x]
	        ye = [ddf_allk[i][name] / P.beta_report for i in x]
	        line = ax.bar(x + idx * width, y, width, color=colors[name], yerr=ye, 
	                      lw=0.1 * elw, error_kw=dict(elinewidth=elw, ecolor='black', capsize=0.5 * elw))
	        lines.append(line[0])

	    ax.set_xlabel('States', fontsize=12, color='#151B54')
	    ax.set_ylabel(f'$\Delta G$ {P.units}', fontsize=12, color='#151B54')
	    ax.set_xticks(x + 0.5 * width * len(P.methods), 
	                  tuple([f'{i}--{i+1}' for i in x]), fontsize=8)
	    ax.set_yticklabels(ax.get_yticks(), fontsize=8)
	    ax.set_xlim([x[0], x[-1] + len(lines) * width])

	    for dir in ['right', 'top', 'bottom']:
	        ax.spines[dir].set_color('none')
	    ax.yaxis.set_ticks_position('left')
	    ax.get_xticklines().set_visible(False)

	    leg = ax.legend(lines, tuple(P.methods), loc=3, ncol=2, prop=FP(size=10), fancybox=True)
	    leg.get_frame().set_alpha(0.5)

	    ax.set_title('The free energy change breakdown', fontsize=12)
	    plt.savefig(os.path.join(P.output_directory, 'dF_state_long.pdf'), bbox_inches='tight')
	    plt.close(fig)

	def plot_dF_vs_lambda2(nb=10):
	    """Plots the free energy differences evaluated for each pair of adjacent states for all methods.
	    The layout is approximately 'nb' bars per subplot."""
	    x = np.arange(len(df_allk))
	    if len(x) < nb:
	        return

	    xs = np.array_split(x, len(x) / nb + 1)
	    mnb = max(len(i) for i in xs)

	    fig, axs = plt.subplots(len(xs), 1, figsize=(8, 6))

	    width = 1. / (len(P.methods) + 1)
	    elw = 30 * width

	    colors = {'TI': '#C45AEC', 'TI-CUBIC': '#33CC33', 'DEXP': '#F87431', 'IEXP': '#FF3030', 'GINS': '#EAC117',
	              'GDEL': '#347235', 'BAR': '#6698FF', 'UBAR': '#817339', 'RBAR': '#C11B17', 'MBAR': '#F9B7FF'}

	    for ndx, x in enumerate(xs, start=1):
	        lines = []
	        ax = axs[ndx - 1]
	        for idx, name in enumerate(P.methods):
	            y = [df_allk[i][name] / P.beta_report for i in x]
	            ye = [ddf_allk[i][name] / P.beta_report for i in x]
	            line = ax.bar(x + idx * width, y, width, color=colors[name], yerr=ye,
	                          lw=0.05 * elw, error_kw=dict(elinewidth=elw, ecolor='black', capsize=0.5 * elw))
	            lines.append(line[0])

	        for dir in ['left', 'right', 'top', 'bottom']:
	            if dir == 'left':
	                ax.yaxis.set_ticks_position(dir)
	            else:
	                ax.spines[dir].set_color('none')

	        ax.set_yticklabels(ax.get_yticks(), fontsize=10)
	        ax.xaxis.set_ticks([])
	        for i in x + 0.5 * width * len(P.methods):
	            ax.annotate(f'$\mathrm{{{i}-{i+1}}}$', xy=(i, 0), xycoords=('data', 'axes fraction'),
	                        xytext=(0, -2), size=10, textcoords='offset points', va='top', ha='center')
	        ax.set_xlim([x[0], x[-1] + len(lines) * width + (mnb - len(x))])

	    leg = ax.legend(lines, tuple(P.methods), loc=0, ncol=2, prop=FP(size=8),
	                    title=f'$\mathrm{{\Delta G\/{P.units}\/}}\mathit{{vs.}}\/\mathrm{{lambda\/pair}}$', fancybox=True)
	    leg.get_frame().set_alpha(0.5)

	    plt.savefig(os.path.join(P.output_directory, 'dF_state.pdf'), bbox_inches='tight')
	    plt.close(fig)

	def plot_TI():
	    """Plots the ave_dhdl array as a function of the lambda value.
	    If (TI and TI-CUBIC in methods) -- plots the TI integration area and the TI-CUBIC interpolation curve,
	    elif (only one of them in methods) -- plots the integration area of the method."""
	    min_dl = dlam[dlam != 0].min()
	    S = int(0.4 / min_dl)

	    fig, ax = plt.subplots(figsize=(8, 6))

	    ax.spines['bottom'].set_position('zero')
	    ax.spines['top'].set_color('none')
	    ax.spines['right'].set_color('none')
	    ax.xaxis.set_ticks_position('bottom')
	    ax.yaxis.set_ticks_position('left')

	    for k, spine in ax.spines.items():
	        spine.set_zorder(12.2)

	    xs, ndx, dx = [0], 0, 0.001
	    colors = ['r', 'g', '#7F38EC', '#9F000F', 'b', 'y']
	    min_y, max_y = 0, 0

	    lines = []
	    lv_names2 = [r'$%s$' % P.lv_names[j].capitalize() for j in range(n_components) if not (ave_dhdl[:, j] == 0).all()]

	    for j in range(n_components):

	        y = ave_dhdl[:, j]
	        if not (y == 0).all():

	            # Get the coordinates.
	            lj = lchange[:, j]
	            x = lv[:, j][lj]
	            y = y[lj] / P.beta_report

	            if 'TI' in P.methods:
	                # Plot the TI integration area.
	                ss = 'TI'
	                for i in range(len(x) - 1):
	                    min_y = min(y.min(), min_y)
	                    max_y = max(y.max(), max_y)
	                    if i % 2 == 0:
	                        ax.fill_between(x[i:i + 2] + ndx, 0, y[i:i + 2], color=colors[ndx], alpha=1.0)
	                    else:
	                        ax.fill_between(x[i:i + 2] + ndx, 0, y[i:i + 2], color=colors[ndx], alpha=0.5)

	                ax.plot([-100 * wnum for wnum in range(len(lv_names2))], 
	                        [0 * wnum for wnum in range(len(lv_names2))], 
	                        ls='-', color=colors[ndx], label=lv_names2[ndx])

	                if 'TI-CUBIC' in P.methods and not cubspl[j] == 0:
	                    # Plot the TI-CUBIC interpolation curve.
	                    ss += ' and TI-CUBIC'
	                    xnew = np.arange(0, 1 + dx, dx)
	                    ynew = cubspl[j].interpolate(y, xnew)
	                    min_y = min(ynew.min(), min_y)
	                    max_y = max(ynew.max(), max_y)
	                    ax.plot(xnew + ndx, ynew, color='#B6B6B4', ls='-', solid_capstyle='round', lw=3.0)

	            else:
	                # Plot the TI-CUBIC integration area.
	                ss = 'TI-CUBIC'
	                for i in range(len(x) - 1):
	                    xnew = np.arange(x[i], x[i + 1] + dx, dx)
	                                        ynew = cubspl[j].interpolate(y, xnew)
                    ynew[0], ynew[-1] = y[i], y[i+1]
                    min_y = min(ynew.min(), min_y)
                    max_y = max(ynew.max(), max_y)
                    if i % 2 == 0:
                        ax.fill_between(xnew + ndx, 0, ynew, color=colors[ndx], alpha=1.0)
                    else:
                        ax.fill_between(xnew + ndx, 0, ynew, color=colors[ndx], alpha=0.5)

	            # Store the abscissa values and update the subplot index.
	            xs += (x + ndx).tolist()[1:]
	            ndx += 1

	    # Make sure the tick labels are not overcrowded.
	    xs = np.array(xs)
	    dl_mat = np.array([xs - i for i in xs])
	    ri = range(len(xs))

	    def getInd(r=ri, z=[0]):
	        primo = r[0]
	        min_dl = ndx * 0.02 * 2 ** (primo > 10)
	        if dl_mat[primo].max() < min_dl:
	            return z
	        for i in r:
	            for j in range(len(xs)):
	                if dl_mat[i, j] > min_dl:
	                    z.append(j)
	                    return getInd(ri[j:], z)

	    xt = [i if (i in getInd()) else '' for i in range(K)]
	    plt.xticks(xs[1:], xt[1:], fontsize=10)
	    plt.yticks(fontsize=10)

	    for tick in ax.get_xticklines():
	        tick.set_visible(False)
	    plt.xlim(0, ndx)
	    min_y *= 1.01
	    max_y *= 1.01
	    plt.ylim(min_y, max_y)

	    for i, j in zip(xs[1:], xt[1:]):
	        plt.annotate(('%.2f' % (i - 1.0 if i > 1.0 else i) if not j == '' else ''), xy=(i, 0), xytext=(i, 0.01), size=10, rotation=90, textcoords=('data', 'axes fraction'), va='bottom', ha='center', color='#151B54')
	    if ndx > 1:
	        lenticks = len(ax.get_ymajorticklabels()) - 1
	        if min_y < 0: lenticks -= 1
	        if lenticks < 5:
	            from matplotlib.ticker import AutoMinorLocator as AML
	            ax.yaxis.set_minor_locator(AML())
	    plt.grid(which='both', color='w', lw=0.25, axis='y', zorder=12)
	    plt.ylabel(r'$\mathrm{\langle{\frac{ \partial U } { \partial \lambda }}\rangle_{\lambda}\/%s}$' % P.units, fontsize=20, color='#151B54')
	    plt.annotate('$\mathit{\lambda}$', xy=(0, 0), xytext=(0.5, -0.05), size=18, textcoords='axes fraction', va='top', ha='center', color='#151B54')
	    if not P.software.title() == 'Sire':
	        lege = ax.legend(prop=FP(size=14), frameon=False, loc=1)
	        for l in lege.legendHandles:
	            l.set_linewidth(10)
	    plt.savefig(os.path.join(P.output_directory, 'dhdl_TI.pdf'))
	    plt.close(fig)
    	return

    print("Plotting the free energy breakdown figure...")
	plotdFvsLambda1()
	plotdFvsLambda2()
	if ('TI' in P.methods or 'TI-CUBIC' in P.methods):
	    print("Plotting the TI figure...")
	    plotTI()

#===================================================================================================
# MAIN
#===================================================================================================

def main():
    global dhdlt
    global u_klt
    global simulation_profile
    global K
    global n_components
    global pymbar
    global dhdl
    global N_k
    global lv
    global dlam
    global ave_dhdl
    global std_dhdl
    global lchange
    global cubspl
    global mapl
    global u_kln
    global Deltaf_ij
    global dDeltaf_ij
    global df_allk
    global ddf_allk
    global nsnapshots
    global pl
    global FP
    global matplotlib

    # Timing.
    start_time = time.time()
    print("Started on", time.asctime())
    print('Command line was:', ' '.join(sys.argv), '\n')

    # Simulation profile to gather information about the simulation.
    simulation_profile = parser.parse_args()[0]

    simulation_profile.methods = ['TI','TI-CUBIC','DEXP','IEXP','BAR','MBAR']
    simulation_profile.units, simulation_profile.beta, simulation_profile.beta_report = check_units_and_more(simulation_profile.units)

    if (numpy.array([simulation_profile.bForwrev, simulation_profile.breakdown, simulation_profile.bCFM, simulation_profile.overlap]) != 0).any():
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.font_manager import FontProperties as FP

    if simulation_profile.software.title() == 'Gromacs':
        import gromacs_parser
        nsnapshots, lv, dhdlt, u_klt = gromacs_parser.read_data_gromacs(simulation_profile)
    else:
        lineno = currentframe().f_lineno
        print("\n\n%s\n Looks like there is no proper parser available to process your files. \n Please modify lines %d and %d of this script.\n%s\n\n" % (78 * "*", lineno + 3, lineno + 4, 78 * "*"))
        #### LINES TO BE MODIFIED
        # import YOUR_OWN_FILE_PARSER
        # nsnapshots, lv, dhdlt, u_klt = YOUR_OWN_FILE_PARSER.your_data_parser(*args, **kwargs)
        #### All the four are numpy arrays.
        #### lv           is the array of lambda vectors.
        #### nsnapshots   is the number of equilibrated snapshots per each state.
        #### dhdlt[k,n,t] is the derivative of energy component n with respect to state k of snapshot t
        #### u_klt[k,m,t] is the reduced potential energy of snapshot t of state k evaluated at state m

    K, n_components = lv.shape
    lchange = get_lchange(lv)

    # Check for all zeros in data files
    all_zeros = not numpy.any(dhdlt) and not numpy.any(u_klt)

    if all_zeros:
        print("WARNING: Found all zeros in input data.")
        zero_output(K, simulation_profile)
        if simulation_profile.bForwrev:
            zero_dFt(K, simulation_profile, nsnapshots)
        sys.exit("Exiting...")

    if (numpy.array(['Sire', 'Gromacs', 'Amber']) == simulation_profile.software.title()).any():
        dhdl, N_k, u_kln = uncorrelate(sta=numpy.zeros(K, int), fin=nsnapshots, do_dhdl=True)
    elif simulation_profile.software.title() == 'Desmond':
        N_k, u_kln = uncorrelate(sta=numpy.zeros(K, int), fin=nsnapshots, do_dhdl=False)

    # Estimate free energy difference with MBAR -- all states at once.
    if 'MBAR' in simulation_profile.methods:
        Deltaf_ij, dDeltaf_ij = estimate_with_MBAR(u_kln, N_k, simulation_profile.relative_tolerance, regular_estimate=True)

    # The TI preliminaries.
    if 'TI' in simulation_profile.methods or 'TI-CUBIC' in simulation_profile.methods:
        dlam, ave_dhdl, std_dhdl = TI_prelim(lv)
    if 'TI-CUBIC' in simulation_profile.methods:
        cubspl, mapl = get_splines(lchange)

    # Call other methods, print stats, store results.
    df_allk, ddf_allk = estimate_pairs()
    total_energies()

    # Plot figures.
    if simulation_profile.bForwrev:
        dF_t()
    if simulation_profile.breakdown:
        plot_dFvsLambda()

    print("\nTime spent: %s hours, %s minutes, and %s seconds.\nFinished on %s" % time_statistics(start_time))

if __name__ == "__main__":
	main()

#===================================================================================================
#                                   End of the script
#===================================================================================================

