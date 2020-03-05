#!/usr/bin/env python3


#------------------------------------------------------------------------------------
# Compute the L1 norm for all given values of dt and plot it.
# Also compute slope of the plot points for all given methods/limiters used.
# You need to define the limiters and dts used below.
#
# You also need to define which "shape" the output filenames have,
# and how many cells there were, as the first two cmdline args.
#------------------------------------------------------------------------------------


# first things first: check whether you can import the hydro python modules
from check_module_is_in_pythonpath import try_to_import
try_to_import()

import numpy as np
import matplotlib
from matplotlib import pyplot as plt

from hydro_utils import get_only_cmdlinearg, get_all_files_with_same_basename
from hydro_io import read_output, read_ic, check_file_exists

from sys import argv
from scipy import stats


shape = argv[1]
nxstr = argv[2]
nx = int(nxstr)


# Plot parameters
params = {
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'font.size': 12,
    'font.family': 'serif',
    'font.serif': 'DejaVu Sans',
    'legend.fontsize': 8,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'figure.dpi' : 200,
    'lines.markersize' : 6,
    'lines.linewidth' : 1.,
    'figure.subplot.left'    : 0.16,
    'figure.subplot.right'   : 0.97,
    'figure.subplot.bottom'  : 0.12,
    'figure.subplot.top'     : 0.92,
}

matplotlib.rcParams.update(params)




def L1norm(expected, obtained):
    """
    Compute the L1 norm for obtained results by comparing it with the
    expected results
    exptected, obtained are both supposed to be np arrays of equal shape,
    here written for 1D
    """
    nx = expected.shape[0]
    L1 = np.sum(np.absolute(expected - obtained)) / nx
    return L1




dt = [0.0002, 0.0001, 0.00009, 0.00007, 0.00005, 0.00003, 0.00001]

limiters = ['NO_LIMITER', 'MINMOD', 'MC', "VANLEER", "SUPERBEE"]
limiter_names = ["no limiter", "minmod", 'MC', "van Leer", "superbee"]

file_prefixes = ['ADVECTION_PWCONST-'] + ['ADVECTION_PWLIN-'+limiter+"-" for limiter in limiters]
labels = ['pcw const'] + limiter_names



dx = 1./nx
x = np.array([(i+0.5)*dx for i in range(nx)])
dt = np.array(dt)
dtlog = np.log(dt)



def analytical_solution(shape, x):
    """
    compute the analytical solution for a given shape and
    at positions x

    must be the same as in your IC, obviously
    """

    for i in range(nx):
        while x[i] < 0:
            x[i] += 1
        while x[i] > 1:
            x[i] -= 1

    sol = np.ones(nx, dtype=np.float)

    if shape == "gaussian":
        for i in range(nx):
            sol[i] = 1 + 2 * np.exp(-((x[i]-0.5)**2/0.05))
    elif shape == "step":
        for i in range(nx):
            if x[i] > 1./3 and x[i] < 2./3:
                sol[i] = 2

    return sol





if __name__ == "__main__":

    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111)
   

    for f in range(len(file_prefixes)):

        accuracy = np.zeros(len(dt))

        for i,c in enumerate(dt):
            fname = file_prefixes[f]+"{0:.8f}-".format(c)+shape+"-"+str(nx)+"-0001.out"

            ndim, rho, u, p, t, nstep = read_output(fname)

            xsol = x - u*t
            expected = analytical_solution(shape, xsol)

            accuracy[i] = L1norm(expected, rho)

        accuracylog = np.log(np.array(accuracy))
        slope, intercept, r_value, p_value, std_err = stats.linregress(dtlog, accuracylog)

        xp = np.array([dt[0], dt[-1]])
        y = np.exp(slope * np.log(xp) + intercept)
        ax.loglog(xp, y,
            c='C'+str(f), label=labels[f]+", slope = {0:6.3f}".format(slope))
        ax.scatter(dt, accuracy,
            c='C'+str(f))



    ax.legend(loc='upper left')
    ax.set_title("nx = {0:d}, nstep = {1:d}".format(nx, nstep))
    ax.set_xlabel(r"$\Delta t$")
    ax.set_ylabel(r"$\frac{1}{N} \sum\limits_{i=1}^{N} |\rho_i - \rho_{i, IC}|$")




    plt.savefig("limiter_accuracy_dt-"+str(shape)+"-"+str(nx)+".png")
