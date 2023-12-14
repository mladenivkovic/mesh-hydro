#!/usr/bin/env python3

# ------------------------------------------------------------------------------------
# Meant to be called by the /sh/compare-slope-limiters-on-1D-advection.sh script
# ------------------------------------------------------------------------------------


# first things first: check whether you can import the hydro python modules
from check_module_is_in_pythonpath import try_to_import

try_to_import()

import numpy as np
import matplotlib
from matplotlib import pyplot as plt

from mesh_hydro_utils import get_only_cmdlinearg, get_all_files_with_same_basename
from mesh_hydro_io import read_output, read_ic


# Plot parameters
params = {
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "font.size": 12,
    "font.family": "serif",
    "font.serif": "DejaVu Sans",
    "legend.fontsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "text.usetex": True,
    "figure.subplot.left": 0.05,
    "figure.subplot.right": 0.97,
    "figure.subplot.bottom": 0.12,
    "figure.subplot.top": 0.92,
    "figure.subplot.wspace": 0.25,
    "figure.subplot.hspace": 0.25,
    "figure.dpi": 200,
    "lines.markersize": 6,
    "lines.linewidth": 1.0,
}

matplotlib.rcParams.update(params)


def plot_1D_density_only(rho, t=0, label=None, dashed=False, fig=None):
    """
    Create a plot from 1D data. Only plot density.

    rho, u, p:  np arrays of physical quantities
    t:          time of the simulation
    label:      label for this line that we're drawing
    dashed:     whether to draw this line as a black dashed line
    fig:        a pyplot.figure object. If present, plots will be added to the axes of the figure.
                if not, a new one will be generated and returned.

    returns:
        fig:     pyplot.figure() object containing the plots
    """

    if fig is None:
        fig = plt.figure(figsize=(6, 5))
        ax1 = fig.add_subplot(1, 1, 1)
    else:
        axes = fig.axes
        ax1 = axes[0]

    nx = rho.shape[0]
    x = np.linspace(0, 1, nx)

    if dashed:
        ax1.plot(x, rho, "--", label=label, c="k")
    else:
        ax1.plot(x, rho, label=label)
    ax1.set_ylabel("density")

    ax1.set_xlabel("x")
    ax1.set_xlim(0, 1)
    ax1.legend(prop={"size": 6})

    plt.subplots_adjust(left=0.15)
    ax1.set_title(r"$t = {0:7.3f}$, nx $= {1:d}$".format(t, nx))

    return fig


limiters = ["NO_LIMITER", "MINMOD", "MC", "VAN_LEER", "SUPERBEE"]
limiter_names = ["no limiter", "minmod", "MC", "van Leer", "superbee"]

if __name__ == "__main__":
    filelist = []
    for limiter in limiters:
        filelist.append(
            "./advection-1D-four-shapes-ADVECTION_PWLIN-" + limiter + "-1D-0001.out"
        )

    fig = None
    for i, f in enumerate(filelist):
        ndim, rho, u, p, t, nsteps = read_output(f)

        if ndim != 1:
            print("I can't overplot 2D stuff...")
            quit(1)
        else:
            fig = plot_1D_density_only(
                rho, t=t, label=limiter_names[i], dashed=False, fig=fig
            )

    # plot ICs
    ndim, twostate, rho, u, p = read_ic("./advection-1D-four-shapes.dat")
    fig = plot_1D_density_only(
        rho, t=t, label="analytical solution", dashed=True, fig=fig
    )

    plt.savefig("limiter_comparison.png")
