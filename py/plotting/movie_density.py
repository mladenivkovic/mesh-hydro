#!/usr/bin/env python3


# ------------------------------------------------------------------------------------
# Create a plot of all output files with the same base name as given file in
# this directory. Basename in this case means everything before _XXXX.out
# If more than one cmdline arg is given, it will instead interpret every arg as an
# individual file, and only plot those.
# You might want to tweak y axis limits for 1D or colorbar limits in 2D manually at
# the beginning of this script.
#
# Usage:
#   movie_density.py file.out     # will find files with same basename
# or:
#   movie_denisty.py <file1> <file2> ... <file N>
# ------------------------------------------------------------------------------------


# first things first: check whether you can import the hydro python modules
from check_module_is_in_pythonpath import try_to_import

try_to_import()


from mesh_hydro_utils import get_only_cmdlinearg, get_all_files_with_same_basename
from mesh_hydro_io import read_output
from mesh_hydro_plotting import get_figname
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

from sys import argv

import numpy as np
from matplotlib import pyplot as plt


dots = False
file_format = "png"

ymin = 0.9
ymax = 2.1
cbarmin = 0.9
cbarmax = 2.5


def plot_1D_density_only(
    rho, fname, dots=False, t=0, draw_legend=False, nosave=False, fig=None
):
    """
    Create a plot from 1D data. Only plot density.

    rho:         np arrays of physical quantities
    fname:       filename of the data you are plotting. Will be used to generate image filename
    dots:        whether to overplot points over lines.
    t:           time of the simulation
    draw_legend: whether to draw a legend. The line labels will be the time.
    nosave:      don't save this figure.
    fig:         a pyplot.figure object. If present, plots will be added to the axes of the figure.
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

    ax1.plot(x, rho, label="t = {0:7.3f}".format(t))
    if dots:
        ax1.scatter(x, rho)
    ax1.set_ylabel("density")

    ax1.set_xlabel("x")
    ax1.set_xlim(0, 1)
    ax1.set_ylim(ymin, ymax)

    if draw_legend:
        ax1.legend(prop={"size": 6})

    plt.figtext(0.5, 0.95, "t = {0:7.3f}".format(t))

    if not nosave:
        plt.subplots_adjust(left=0.10)  # ATTENTION: This is different in the module

        # if we have more than one line, this figure is used to overplot things.
        # in that case, give it a different name.
        case = "density"
        threshold = 1
        if dots:
            threshold += 1
        if len(ax1.get_lines()) > threshold:
            case = "density"
        figname = get_figname(fname, case=case)

        plt.savefig(figname, format=file_format)
        print("Saved figure", figname)

    plt.close()

    return fig


def plot_2D_density_only(rho, fname, t=0):
    """
    Create a plot from 2D data. Plots density only.

    rho:         np arrays of physical quantities
    fname:       filename of the data you are plotting. Will be used to generate image filename
    """

    fig = plt.figure(figsize=(6, 5))

    nx = rho.shape[0]

    ax1 = fig.add_subplot(1, 1, 1)
    im1 = ax1.imshow(
        rho, origin="lower", extent=(0, 1, 0, 1), vmin=cbarmin, vmax=cbarmax
    )
    ax1.set_ylabel("density")

    ax1.set_xlabel("x")
    ax1.set_xlim(0, 1)
    ax1.set_ylabel("y")
    ax1.set_ylim(0, 1)

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im1, cax=cax)

    plt.figtext(0.5, 0.95, "t = {0:7.3f}".format(t))

    figname = get_figname(fname, case="density")

    plt.savefig(figname, format=file_format)
    print("Saved figure", figname)

    plt.close()

    return


if __name__ == "__main__":
    if len(argv) == 2:
        fname = get_only_cmdlinearg()
        filelist = get_all_files_with_same_basename(fname)
    else:
        filelist = argv[1:]

    for fname in filelist:
        ndim, rho, u, p, t, step = read_output(fname)

        if ndim == 1:
            plot_1D_density_only(rho, fname, dots=False, t=t)
        elif ndim == 2:
            plot_2D_density_only(rho, fname, t=t)
