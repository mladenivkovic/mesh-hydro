#!/usr/bin/env python3


#------------------------------------------------------------------------------------
# Compute the L1 norm for all given values of Ccfl and plot it.
# Also compute slope of the plot points for all given methods used.
# You need to define riemann solvers and Ccfl values manually below.
#
# You need to give the initial conditions file used as the first cmdlinearg.
#------------------------------------------------------------------------------------


# first things first: check whether you can import the hydro python modules
from check_module_is_in_pythonpath import try_to_import
try_to_import()

import numpy as np
import matplotlib
from matplotlib import pyplot as plt

from mesh_hydro_utils import get_only_cmdlinearg, get_all_files_with_same_basename
from mesh_hydro_io import read_output, read_ic, check_file_exists
from mesh_hydro_riemann import riemann_solver

from sys import argv
from scipy import stats

icfile = argv[1]


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




ccfl = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]

solvers = ['EXACT', 'HLLC', 'TRRS', "TSRS"]

file_prefixes = ['GODUNOV-'+ s for s in solvers]
labels = solvers


ccfllog = np.log(np.array(ccfl))




if __name__ == "__main__":

    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111)
   
    for f in range(len(file_prefixes)):

        accuracy = np.zeros(len(ccfl))

        for i,c in enumerate(ccfl):
            fname = file_prefixes[f]+"-{0:.2f}-0001.out".format(c)

            ndim, rho, u, p, t, nstep = read_output(fname)
            nx = rho.shape[0]
            
            _, _, rhoIC, uIC, pIC = read_ic(icfile, nx)
            rho_sol, u_sol, p_sol = riemann_solver(rhoIC, uIC, pIC, t)

            accuracy[i] = L1norm(rho_sol, rho)

        ax.loglog(ccfl, accuracy, 
            c='C'+str(f), label=labels[f])
        ax.scatter(ccfl, accuracy, 
            c='C'+str(f))

    # add a line with slope 0.5 starting at the lowest point for comparison
    xp1 = ccfl[-1]
    yp1 = accuracy[-1]
    xp2 = ccfl[0]
    slope = 0.5
    yp2 = np.exp((slope * np.log(xp2/xp1)) + np.log(yp1))

    ax.plot([xp1, xp2], [yp1, yp2], 'k:', label="slope 0.5")


    ax.legend(loc='upper left')
    ax.set_title("nx = {0:d}, nsteps = {1:d}".format(nx, nstep))
    ax.set_xlabel(r"$C_{cfl}$")
    ax.set_ylabel(r"$\frac{1}{N} \sum\limits_{i=1}^{N} |\rho_i - \rho_{i, IC}|$")

    plt.savefig("godunov-convergence-CFL-"+str(nx)+".png")
