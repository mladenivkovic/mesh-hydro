#!/usr/bin/env python3


# ------------------------------------------------------------------------------------
# Compute the L1 norm for all given values of nx and plot it.
# Also compute slope of the plot points for all given methods/limiters used.
# You need to define the Riemann solvers and nx values used below.
#
# You need to give the initial conditions file used as the first cmdlinearg.
# ------------------------------------------------------------------------------------

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
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "font.size": 12,
    "font.family": "serif",
    "font.serif": "DejaVu Sans",
    "legend.fontsize": 8,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "text.usetex": True,
    "figure.dpi": 200,
    "lines.markersize": 6,
    "lines.linewidth": 1.0,
    "figure.subplot.left": 0.16,
    "figure.subplot.right": 0.97,
    "figure.subplot.bottom": 0.12,
    "figure.subplot.top": 0.92,
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


nx = [5000, 10000, 20000, 50000, 80000, 100000]
dx = [1 / n for n in nx]

solvers = ["EXACT", "HLLC", "TRRS", "TSRS"]
file_prefixes = ["GODUNOV-" + s for s in solvers]
labels = solvers

dxlog = np.log(dx)


if __name__ == "__main__":
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    for f in range(len(file_prefixes)):
        accuracy = np.zeros(len(dx))

        for i, c in enumerate(dx):
            fname = file_prefixes[f] + "-" + str(nx[i]) + "-0001.out".format(c)

            ndim, rho, u, p, t, nstep = read_output(fname)

            _, _, rhoIC, uIC, pIC = read_ic(icfile, nx[i])
            rho_sol, u_sol, p_sol = riemann_solver(rhoIC, uIC, pIC, t)

            accuracy[i] = L1norm(rho_sol, rho)

        accuracylog = np.log(np.array(accuracy))
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            dxlog, accuracylog
        )

        xp = np.array([dx[0], dx[-1]])
        yp = np.exp(slope * np.log(xp) + intercept)
        ax.loglog(
            xp, yp, c="C" + str(f), label=labels[f] + ", slope = {0:6.3f}".format(slope)
        )
        ax.scatter(dx, accuracy, c="C" + str(f))

    ax.legend(loc="upper left")
    ax.set_title("nsteps = {1:d}".format(nx, nstep))
    ax.set_xlabel(r"$\Delta x$")
    ax.set_ylabel(r"$\frac{1}{N} \sum\limits_{i=1}^{N} |\rho_i - \rho_{i, IC}|$")

    plt.savefig("godunov-convergence-dx.png")
