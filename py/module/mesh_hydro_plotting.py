#!/usr/bin/env python3


# ---------------------------
# Plotting functions
# ---------------------------


import matplotlib
from matplotlib import pyplot as plt

#  matplotlib.use("Agg")
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from mpl_toolkits.mplot3d import Axes3D

import numpy as np


# Plot parameters
params = {
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "font.size": 14,
    "font.family": "serif",
    "font.serif": "DejaVu Sans",
    "legend.fontsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "text.usetex": True,
    "figure.subplot.left": 0.05,
    "figure.subplot.right": 0.97,
    "figure.subplot.bottom": 0.12,
    "figure.subplot.top": 0.92,
    "figure.subplot.wspace": 0.25,
    "figure.subplot.hspace": 0.25,
    "figure.dpi": 200,
    "lines.markersize": 6,
    "lines.linewidth": 2.0,
}

matplotlib.rcParams.update(params)

# set file format here if necessary. Allowed options: png, pdf
file_format = "png"
#  file_format = "pdf"


def plot_1D(rho, u, p, draw_legend=False, fig=None, kwargs={}):
    """
    Create a plot from 1D data.

    rho, u, p:   np arrays of physical quantities
    draw_legend: whether to draw a legend.
    fig:         a pyplot.figure object. If present, plots will be added to the axes of the figure.
                 if not, a new one will be generated and returned.

    kwargs get passed to matplotlib.pyplot.plot(), and need to be a dictionnary

    returns:
        fig:     pyplot.figure() object containing the plots
    """

    if fig is None:
        fig = plt.figure(figsize=(15, 5))
        ax1 = fig.add_subplot(1, 3, 1)
        ax2 = fig.add_subplot(1, 3, 2)
        ax3 = fig.add_subplot(1, 3, 3)
    else:
        axes = fig.axes
        ax1 = axes[0]
        ax2 = axes[1]
        ax3 = axes[2]

    nx = rho.shape[0]
    x = np.linspace(0, 1, nx)

    ax1.plot(
        x, rho, **kwargs,
    )
    ax1.set_ylabel("density")

    ax2.plot(
        x, u, **kwargs,
    )
    ax2.set_ylabel("velocity")

    ax3.plot(
        x, p, **kwargs,
    )
    ax3.set_ylabel("pressure")

    for ax in fig.axes:
        ax.set_xlabel("x")
        ax.set_xlim(0, 1)
        if draw_legend:
            ax.legend(prop={"size": 6})

    return fig


def plot_1D_density_only(rho, draw_legend=False, fig=None, kwargs={}):
    """
    Create a plot from 1D data. Only plot density.

    rho, u, p:   np arrays of physical quantities
    draw_legend: whether to draw a legend.
    fig:         a pyplot.figure object. If present, plots will be added to the axes of the figure.
                 if not, a new one will be generated and returned.

    kwargs get passed to matplotlib.pyplot.plot(), and need to be a dictionnary

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

    ax1.plot(
        x, rho, **kwargs,
    )
    ax1.set_ylabel("density")

    ax1.set_xlabel("x")
    ax1.set_xlim(0, 1)

    if draw_legend:
        ax1.legend(prop={"size": 6})

    plt.subplots_adjust(left=0.15)

    return fig


def plot_2D_density_only(rho, t=None, kwargs={}):
    """
    Create a plot from 2D data. Plots density only.

    rho:         np arrays of physical quantities
    fname:       filename of the data you are plotting. Will be used to generate image filename
    t:           time of simulation, optional. Will be put on the plot to label it. If it is a string, it will be used just as the label.

    kwargs get passed to matplotlib.pyplot.imshow(), and need to be a dictionnary

    returns:
        fig:     pyplot.figure() object containing the plots
    """

    fig = plt.figure(figsize=(6, 5))

    nx = rho.shape[0]

    ax1 = fig.add_subplot(1, 1, 1)
    im1 = ax1.imshow(rho, origin="lower", extent=(0, 1, 0, 1), **kwargs,)
    ax1.set_title("density")

    ax1.set_xlabel("x")
    ax1.set_xlim(0, 1)
    ax1.set_ylabel("y")
    ax1.set_ylim(0, 1)

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im1, cax=cax)

    if t is not None:
        if isinstance(t, str):
            text = t
        elif isinstance(t, float):
            text = r"$t = ${0:.3f}".format(t)
        else:
            raise ValueError(
                "Got weird data type for label (t). t=", t, "type(t)=", type(t)
            )
        plt.figtext(0.05, 0.95, text)

    return fig


def plot_2D(rho, u, p, t=None, kwargs={}):
    """
    Create a plot from 2D data.
    rho, u, p:   np arrays of physical quantities
    t:           time of the simulation
    fig:         a pyplot.figure object. If present, plots will be added to the axes of the figure.
                 if not, a new one will be generated and returned.

    kwargs get passed to matplotlib.pyplot.imshow(), and need to be a dictionnary

    returns:
        fig:     pyplot.figure() object containing the plots
    """

    fig = plt.figure(figsize=(21, 5))

    nx = rho.shape[0]

    ax1 = fig.add_subplot(1, 4, 1)
    im1 = ax1.imshow(rho, origin="lower", extent=(0, 1, 0, 1), **kwargs,)
    ax1.set_title("density")

    ax2 = fig.add_subplot(1, 4, 2)
    im2 = ax2.imshow(u[:, :, 0], origin="lower", extent=(0, 1, 0, 1), **kwargs,)
    ax2.set_title("velocity in x direction")

    ax3 = fig.add_subplot(1, 4, 3)
    im3 = ax3.imshow(u[:, :, 1], origin="lower", extent=(0, 1, 0, 1), **kwargs,)
    ax3.set_title("velocity in y direction")

    ax4 = fig.add_subplot(1, 4, 4)
    im4 = ax4.imshow(p, origin="lower", extent=(0, 1, 0, 1), **kwargs,)
    ax4.set_title("pressure")

    for ax in fig.axes:
        ax.set_xlabel("x")
        ax.set_xlim(0, 1)
        ax.set_ylabel("y")
        ax.set_ylim(0, 1)

    for im, ax in [(im1, ax1), (im2, ax2), (im3, ax3), (im4, ax4)]:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)

    if t is not None:
        if isinstance(t, str):
            text = t
        elif isinstance(t, float):
            text = r"$t = ${0:.3f}".format(t)
        else:
            raise ValueError(
                "Got weird data type for label (t). t=", t, "type(t)=", type(t)
            )
        plt.figtext(0.05, 0.95, text)

    return fig


def plot_2D_velnorm(rho, u, p, t=None, kwargs={}):
    """
    Create a plot from 2D data. Plot velocity norm, not components.
    rho, u, p:   np arrays of physical quantities
    t:           time of the simulation

    kwargs get passed to matplotlib.pyplot.imshow(), and need to be a dictionnary

    returns:
        fig:    pyplot figure object
    """

    unorm = np.sqrt(u[:, :, 0] ** 2 + u[:, :, 1] ** 2)

    fig = plt.figure(figsize=(16, 5))

    nx = rho.shape[0]

    ax1 = fig.add_subplot(1, 3, 1)
    im1 = ax1.imshow(rho, origin="lower", extent=(0, 1, 0, 1), **kwargs,)
    ax1.set_title("density")

    ax2 = fig.add_subplot(1, 3, 2)
    im2 = ax2.imshow(u[:, :, 0], origin="lower", extent=(0, 1, 0, 1), **kwargs,)
    ax2.set_title("velocity norm")

    ax3 = fig.add_subplot(1, 3, 3)
    im3 = ax3.imshow(p, origin="lower", extent=(0, 1, 0, 1), **kwargs,)
    ax3.set_title("pressure")

    for ax in fig.axes:
        ax.set_xlabel("x")
        ax.set_xlim(0, 1)
        ax.set_ylabel("y")
        ax.set_ylim(0, 1)

    for im, ax in [(im1, ax1), (im2, ax2), (im3, ax3)]:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)

    if t is not None:
        if isinstance(t, str):
            text = t
        elif isinstance(t, float):
            text = r"$t = ${0:.3f}".format(t)
        else:
            raise ValueError(
                "Got weird data type for label (t). t=", t, "type(t)=", type(t)
            )
        plt.figtext(0.05, 0.95, text)

    return fig


def plot_2D_in_3D(rho, u, p, t=None, kwargs={}):
    """
    Create a 3D plot from 2D data.
    rho, u, p:   np arrays of physical quantities
    t:           time of the simulation

    kwargs get passed to matplotlib.pyplot.plot_surface(), and need to be a dictionnary

    returns:
        fig:    pyplot figure object
    """

    fig = plt.figure(figsize=(21, 5))

    nx = rho.shape[0]

    X = np.linspace(0, 1, nx)
    Y = np.linspace(0, 1, nx)
    X, Y = np.meshgrid(X, Y)

    ax1 = fig.add_subplot(1, 4, 1, projection="3d")
    im1 = ax1.plot_surface(X, Y, rho, cmap="viridis", **kwargs,)
    ax1.set_title("density")

    ax2 = fig.add_subplot(1, 4, 2, projection="3d")
    im2 = ax2.plot_surface(X, Y, u[:, :, 0], cmap="viridis", **kwargs,)
    ax2.set_title("velocity in x direction")

    ax3 = fig.add_subplot(1, 4, 3, projection="3d")
    im3 = ax3.plot_surface(X, Y, u[:, :, 1], cmap="viridis", **kwargs,)
    ax3.set_title("velocity in y direction")

    ax4 = fig.add_subplot(1, 4, 4, projection="3d")
    im4 = ax4.plot_surface(X, Y, p, cmap="viridis", **kwargs,)
    ax4.set_title("pressure")

    for ax in fig.axes:
        ax.set_xlabel("x")
        ax.set_xlim(0, 1)
        ax.set_ylabel("y")
        ax.set_ylim(0, 1)

    if t is not None:
        if isinstance(t, str):
            text = t
        elif isinstance(t, float):
            text = r"$t = ${0:.3f}".format(t)
        else:
            raise ValueError(
                "Got weird data type for label (t). t=", t, "type(t)=", type(t)
            )
        plt.figtext(0.05, 0.95, text)

    return fig


def save_plot(fig, fname=None, case=None, fname_force=None):
    """
    Save the figure. If fname and case are given, it will first generate
    a descriptive figure name. Otherwise, it will save it as 'hydro-plot.png'
    """

    if fname_force is not None:
        plt.savefig(fname_force, fig=fig)
        print("Saved figure", fname_force)

    else:

        if fname is not None:

            ax = fig.axes[0]
            if len(ax.get_lines()) > 1:
                # if we have more than one line, this figure is used to overplot things.
                # in that case, give it a different name.
                if case is None:
                    case = "overplotted"
                else:
                    # Note: this skips the case "not-overplotted", which is as intended
                    if "overplotted" not in case:
                        case += "-overplotted"
            fname = get_figname(fname, case)
        else:
            fname = "hydro-plot.png"

        plt.savefig(fname, fig=fig)
        print("Saved figure", fname)

    plt.close()
    return


def get_figname(fname, case=None):
    """
    Generate figure name using initial filename fname.
    Remove the file suffix, if present, and add a png.

    if case is not None, it will add something to the file name
    so it will be distinguishable. 
    Accepted cases:
        "density":              for density-only plotting
        "overplotted":          mutliple plots on one axis
        "density-overplotted":  density only with multiple lines per axis   

    returns:
        figname:    figure name string
    """

    # start from last letter, look for a dot to find the suffix
    # if you reach a slash first, stop there
    nchars = len(fname)
    dotindex = None
    for c in range(nchars):
        if fname[nchars - c - 1] == "/":
            break
        if fname[nchars - c - 1] == ".":
            dotindex = nchars - c - 1
            break

    # now extract the actual file basename
    if dotindex is None:
        figname = fname
    else:
        figname = fname[:dotindex]

    if case is not None:
        if case == "density":
            figname += "-density-only"
        if case == "overplotted":
            # first remove snapshot number
            figname = figname[:-5]
            figname += "-overplotted"
        if case == "density-overplotted":
            figname = figname[:-5]
            figname += "-density-only-overplotted"
        if case == "3D":
            figname += "-3D"
        if case == "riemann-solver":
            figname += "-riemann-solution"
        if case == "not-overplotted":
            pass
        if case == "artsy":
            figname += "-artsy"

    figname += "." + file_format

    return figname
