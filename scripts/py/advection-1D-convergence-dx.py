#!/usr/bin/env python3


# ------------------------------------------------------------------------------------
# Compute the L1 norm for all given values of nx and plot it.
# Also compute slope of the plot points for all given methods/limiters used.
# You need to define limiters and nx values manually below.
#
# You also need to define which "shape" the output filenames have, i.e. what we're
# advecting (four-shapes, gaussian, step), which is taken as first cmd line arg.
# You also need to give Ccfl as second command line arg.
#
# Meant to be run with /sh/get-dx-convergence-1D-advection.sh
# ------------------------------------------------------------------------------------


# first things first: check whether you can import the hydro python modules
from check_module_is_in_pythonpath import try_to_import

try_to_import()

import numpy as np
import matplotlib
from matplotlib import pyplot as plt

from mesh_hydro_utils import get_only_cmdlinearg, get_all_files_with_same_basename
from mesh_hydro_io import read_output, read_ic

from sys import argv
from scipy import stats


shape = argv[1]

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


limiters = ["NO_LIMITER", "MINMOD", "MC", "VANLEER", "SUPERBEE"]
limiter_names = ["no limiter", "minmod", "MC", "van Leer", "superbee"]
nx = [100, 200, 500, 800, 1000, 1500, 2000]

file_prefixes = ["ADVECTION_PWCONST-"] + [
    "ADVECTION_PWLIN-" + limiter + "-" for limiter in limiters
]
labels = ["pcw const"] + limiter_names


dx = np.array([1.0 / x for x in nx])
dxlog = np.log(dx)


def analytical_solution(shape, x):
    """
    compute the analytical solution for a given shape and
    at positions x

    must be the same as you set in your IC, obviously!
    """

    nx = x.shape[0]

    for i in range(nx):
        while x[i] < 0:
            x[i] += 1
        while x[i] > 1:
            x[i] -= 1

    sol = np.ones(nx, dtype=np.float)

    if shape == "gaussian":
        for i in range(nx):
            sol[i] = 1 + 2 * np.exp(-((x[i] - 0.5) ** 2 / 0.05))

    elif shape == "step":
        for i in range(nx):
            if x[i] > 1.0 / 3 and x[i] < 2.0 / 3:
                sol[i] = 2

    else:
        print("Error: shape", shape, "not recognized")
        quit()

    return sol


if __name__ == "__main__":
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    for f in range(len(file_prefixes)):
        accuracy = np.zeros(len(dx))

        for i, c in enumerate(nx):
            fname = file_prefixes[f] + shape + "-" + str(c) + "-0001.out"

            ndim, rho, u, p, t, nstep = read_output(fname)
            x = np.array([(j + 0.5) * dx[i] for j in range(c)])
            xsol = x - u * t
            expected = analytical_solution(shape, xsol)

            accuracy[i] = L1norm(expected, rho)

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

    ax.legend(loc="lower right")
    ax.set_title(r"nsteps = {:d}".format(nstep))
    ax.set_xlabel(r"dx")
    ax.set_ylabel(r"$\frac{1}{N} \sum\limits_{i=1}^{N} |\rho_i - \rho_{i, IC}|$")

    plt.savefig("limiter_accuracy_dx-" + shape + ".png")
