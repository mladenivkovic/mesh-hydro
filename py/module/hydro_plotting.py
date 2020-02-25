#!/usr/bin/env python3


#---------------------------
# Plotting functions
#---------------------------


import matplotlib
from matplotlib import pyplot as plt
#  matplotlib.use("Agg")
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from mpl_toolkits.mplot3d import Axes3D

import numpy as np


# Plot parameters
params = {
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'font.size': 12,
    'font.family': 'serif',
    'font.serif': 'DejaVu Sans',
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'figure.subplot.left'    : 0.05,
    'figure.subplot.right'   : 0.97,
    'figure.subplot.bottom'  : 0.12,
    'figure.subplot.top'     : 0.92,
    'figure.subplot.wspace'  : 0.25,
    'figure.subplot.hspace'  : 0.25,
    'figure.dpi' : 200,
    'lines.markersize' : 6,
    'lines.linewidth' : 2.
}

matplotlib.rcParams.update(params)

# set file format here if necessary. Allowed options: png, pdf
file_format = "png"
#  file_format = "pdf"







def plot_1D(rho, u, p, fname, dots=False, t = None, draw_legend = False, nosave = False, fig = None):
    """
    Create a plot from 1D data.

    rho, u, p:   np arrays of physical quantities
    fname:       filename of the data you are plotting. Will be used to generate image filename
    dots:        whether to overplot points over lines.
    t:           time of the simulation. Used to label lines. If it is a string, the string will be used as the label instead.
    draw_legend: whether to draw a legend. The line labels will be the time.
    nosave:      don't save this figure.
    fig:         a pyplot.figure object. If present, plots will be added to the axes of the figure.
                 if not, a new one will be generated and returned.

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

    text = None
    if t is not None:
        if isinstance(t, str):
            text = t
        elif isinstance(t, float):
            text = r"$t = ${0:.3f}".format(t)
        else:
            raise ValueError("Got weird data type for label (t). t=", t, "type(t)=", type(t))

    ax1.plot(x, rho,
             label = text,
            )
    if dots: ax1.scatter(x, rho)
    ax1.set_ylabel('density')

    ax2.plot(x, u,
             label = text,
            )
    if dots: ax2.scatter(x, u)
    ax2.set_ylabel('velocity')

    ax3.plot(x, p,
             label = text,
            )
    if dots: ax3.scatter(x, p)
    ax3.set_ylabel('pressure')

    for ax in fig.axes:
        ax.set_xlabel("x")
        ax.set_xlim(0,1)
        if draw_legend:
            ax.legend(prop={'size': 6})


    if not nosave:
        # if we have more than one line, this figure is used to overplot things.
        # in that case, give it a different name.
        case = None
        if len(ax.get_lines()) > 1:
            case = "overplotted"
        figname = get_figname(fname, case = case)

        plt.savefig(figname, format=file_format)
        print("Saved figure", figname)

    return fig






def plot_1D_density_only(rho, fname, dots=False, t = 0, draw_legend = False, nosave = False, fig = None):
    """
    Create a plot from 1D data. Only plot density.

    rho, u, p:   np arrays of physical quantities
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

    ax1.plot(x, rho,
             label = "t = {0:7.3f}".format(t),
            )
    if dots: ax1.scatter(x, rho)
    ax1.set_ylabel('density')

    ax1.set_xlabel("x")
    ax1.set_xlim(0,1)
    if draw_legend:
        ax1.legend(prop={'size': 6})

    if not nosave:
        plt.subplots_adjust(left=0.15)

        # if we have more than one line, this figure is used to overplot things.
        # in that case, give it a different name.
        case = "density"
        threshold = 1
        if dots: threshold += 1
        if len(ax1.get_lines()) > threshold:
            case = "density-overplotted"
        figname = get_figname(fname, case = case)

        plt.savefig(figname, format=file_format)
        print("Saved figure", figname)

    return fig








def plot_2D_density_only(rho, fname, t=None):
    """
    Create a plot from 2D data. Plots density only.

    rho:         np arrays of physical quantities
    fname:       filename of the data you are plotting. Will be used to generate image filename
    t:           time of simulation, optional. Will be put on the plot to label it. If it is a string, it will be used just as the label.
    """

    fig = plt.figure(figsize=(6, 5))

    nx = rho.shape[0]

    ax1 = fig.add_subplot(1, 1, 1)
    im1 = ax1.imshow(rho,
            origin='lower', 
            extent=(0,1,0,1),
            )
    ax1.set_title('density')

    ax1.set_xlabel("x")
    ax1.set_xlim(0,1)
    ax1.set_ylabel("y")
    ax1.set_ylim(0,1)


    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im1, cax=cax)

    if t is not None:
        if isinstance(t, str):
            text = t
        elif isinstance(t, float):
            text = r"$t = ${0:.3f}".format(t)
        else:
            raise ValueError("Got weird data type for label (t). t=", t, "type(t)=", type(t))
        plt.figtext(0.05, 0.95, text)


    figname = get_figname(fname, case="density")

    plt.savefig(figname, format=file_format)
    print("Saved figure", figname)

    return fig









def plot_2D(rho, u, p, fname):
    """
    Create a plot from 2D data.
    rho, u, p:   np arrays of physical quantities
    fname:       filename of the data you are plotting. Will be used to generate image filename
    t:           time of the simulation
    draw_legend: whether to draw a legend. The line labels will be the time.
    nosave:      don't save this figure.
    fig:         a pyplot.figure object. If present, plots will be added to the axes of the figure.
                 if not, a new one will be generated and returned.
    """

    fig = plt.figure(figsize=(21, 5))

    nx = rho.shape[0]

    ax1 = fig.add_subplot(1, 4, 1)
    im1 = ax1.imshow(rho,
            origin='lower', 
            extent=(0,1,0,1),
            )
    ax1.set_title('density')

    ax2 = fig.add_subplot(1, 4, 2)
    im2 = ax2.imshow(u[:,:,0],
            origin='lower', 
            extent=(0,1,0,1),
            )
    ax2.set_title('velocity in x direction')

    ax3 = fig.add_subplot(1, 4, 3)
    im3 = ax3.imshow(u[:,:,1],
            origin='lower', 
            extent=(0,1,0,1),
            )
    ax3.set_title('velocity in y direction')

    ax4 = fig.add_subplot(1, 4, 4)
    im4 = ax4.imshow(p, 
            origin='lower', 
            extent=(0,1,0,1),
            )
    ax4.set_title('pressure')

    for ax in fig.axes:
        ax.set_xlabel("x")
        ax.set_xlim(0,1)
        ax.set_ylabel("y")
        ax.set_ylim(0,1)


    for im, ax in [(im1, ax1), (im2, ax2), (im3, ax3), (im4, ax4)]:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)



    figname = get_figname(fname)

    plt.savefig(figname, format=file_format)
    print("Saved figure", figname)

    return fig







def plot_2D_in_3D(rho, u, p, fname):
    """
    Create a 3D plot from 2D data.
    rho, u, p:   np arrays of physical quantities
    fname:       filename of the data you are plotting. Will be used to generate image filename
    t:           time of the simulation
    draw_legend: whether to draw a legend. The line labels will be the time.
    nosave:      don't save this figure.
    fig:         a pyplot.figure object. If present, plots will be added to the axes of the figure.
                 if not, a new one will be generated and returned.
    """

    fig = plt.figure(figsize=(21, 5))

    nx = rho.shape[0]

    X = np.linspace(0, 1, nx)
    Y = np.linspace(0, 1, nx)
    X, Y = np.meshgrid(X, Y)

    ax1 = fig.add_subplot(1, 4, 1, projection='3d')
    im1 = ax1.plot_surface(X, Y, rho, cmap = 'viridis' )
    ax1.set_title('density')

    ax2 = fig.add_subplot(1, 4, 2, projection='3d')
    im2 = ax2.plot_surface(X, Y, u[:,:,0], cmap = 'viridis' )
    ax2.set_title('velocity in x direction')

    ax3 = fig.add_subplot(1, 4, 3, projection='3d')
    im3 = ax3.plot_surface(X, Y, u[:,:,1], cmap = 'viridis' )
    ax3.set_title('velocity in y direction')

    ax4 = fig.add_subplot(1, 4, 4, projection='3d')
    im4 = ax4.plot_surface(X, Y, p, cmap = 'viridis' )
    ax4.set_title('pressure')

    for ax in fig.axes:
        ax.set_xlabel("x")
        ax.set_xlim(0,1)
        ax.set_ylabel("y")
        ax.set_ylim(0,1)

    figname = get_figname(fname, "3D")

    plt.savefig(figname, format=file_format)
    print("Saved figure", figname)

    return fig







def get_figname(fname, case = None):
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
        if fname[nchars-c-1] == '/':
            break
        if fname[nchars-c-1] == ".":
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
            

    figname += "." + file_format

    return figname
