#!/usr/bin/env python3


# -------------------------------------
# Create a nice looking plot without
# any scientific significance.
# plots density only without borders
# nor axes nor colorbars.
#
# NOT INTENDED FOR GENERAL USE.
# -------------------------------------


# first things first: check whether you can import the hydro python modules
from check_module_is_in_pythonpath import try_to_import

try_to_import()

import matplotlib
from matplotlib import pyplot as plt

#  matplotlib.use("Agg")
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from mpl_toolkits.mplot3d import Axes3D

import numpy as np

from mesh_hydro_plotting import plot_1D_density_only, save_plot
from mesh_hydro_io import read_output
from mesh_hydro_utils import get_only_cmdlinearg, label_to_kwargs

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
    "figure.subplot.left": 0.0,
    "figure.subplot.right": 1.0,
    "figure.subplot.bottom": 0.0,
    "figure.subplot.top": 1.0,
    "figure.subplot.wspace": 0.0,
    "figure.subplot.hspace": 0.0,
    "figure.dpi": 200,
    "lines.markersize": 6,
    "lines.linewidth": 2.0,
}

matplotlib.rcParams.update(params)


def plot_2D_density_only(rho, t=None, kwargs={}):
    """
    Create a plot from 2D data. Plots density only.

    rho:         np arrays of physical quantities
    fname:       filename of the data you are plotting. Will be used to generate image filename
    t:           time of simulation, optional. Will be put on the plot to label it. If it is a string, it will be used just as the label.

    kwargs get passed to matplotlib.pyplot.imshow(), and need to be a dictionnary
    """

    fig = plt.figure(figsize=(10, 5))

    nx = rho.shape[0]

    ax1 = fig.add_subplot(1, 1, 1)
    im1 = ax1.imshow(rho, origin="lower", extent=(0, 1, 0, 0.5), **kwargs)

    # turn off axis
    ax1.set_axis_off()

    # cut off margins
    ax1.xaxis.set_major_locator(plt.NullLocator())
    ax1.yaxis.set_major_locator(plt.NullLocator())

    return fig


if __name__ == "__main__":

    fname = get_only_cmdlinearg()
    ndim, rho, u, p, t, step = read_output(fname)

    if ndim == 2:
        fig = plot_2D_density_only(rho, t=t)

        # Upper half for kelvin helmholtz.
        # remember imshow takex array[y, x]
        #  fig = plot_2D_density_only(rho[128:256, :], t=t)

        save_plot(fig, fname, case="artsy")
