#!/usr/bin/env python3


# ------------------------------------------------------------------------------------
# Compute the L1 norm for all given values of Ccfl and plot it.
# You need to define limiters and Ccfl values manually below.
#
# You also need to define which "shape" the output filenames have, as the first cmdlinearg
# ------------------------------------------------------------------------------------

# first things first: check whether you can import the hydro python modules
from check_module_is_in_pythonpath import try_to_import

try_to_import()


import numpy as np
import matplotlib
from matplotlib import pyplot as plt

from mesh_hydro_utils import get_only_cmdlinearg, get_all_files_with_same_basename
from mesh_hydro_io import read_output, read_ic, check_file_exists

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


def analytical_solution(shape, x):
    """
    compute the analytical solution for a given shape and
    at positions x
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


ccfl = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]

limiters = ["NO_LIMITER", "MINMOD", "MC", "VANLEER", "SUPERBEE"]
limiter_names = ["no limiter", "minmod", "MC", "van Leer", "superbee"]

file_prefixes = ["ADVECTION_PWCONST-"] + [
    "ADVECTION_PWLIN-" + limiter + "-" for limiter in limiters
]
labels = ["pcw const"] + limiter_names


ccfllog = np.log(np.array(ccfl))


if __name__ == "__main__":
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    for f in range(len(file_prefixes)):
        accuracy = np.zeros(len(ccfl))

        for i, c in enumerate(ccfl):
            fname = file_prefixes[f] + shape + "-{0:.8f}-0001.out".format(c)

            ndim, rho, u, p, t, nstep = read_output(fname)
            nx = rho.shape[0]
            dx = 1.0 / nx
            x = np.array([(i + 0.5) * dx for i in range(nx)])
            xsol = x - u * t
            expected = analytical_solution(shape, xsol)

            accuracy[i] = L1norm(expected, rho)

        ax.loglog(ccfl, accuracy, c="C" + str(f), label=labels[f])
        ax.scatter(ccfl, accuracy, c="C" + str(f))

    ax.legend(loc="upper left")
    ax.set_title("nx = {0:d}, nsteps = {1:d}".format(nx, nstep))
    ax.set_xlabel(r"$C_{cfl}$")
    ax.set_ylabel(r"$\frac{1}{N} \sum\limits_{i=1}^{N} |\rho_i - \rho_{i, IC}|$")

    plt.savefig("limiter_accuracy_CFL-" + shape + "-" + str(nx) + ".png")
