#!/usr/bin/env python3


# ------------------------------------------------------------
# Create a plot of an output file.
# If it's a 1D output file, it creates the same as the
# plot_result.py script. If it's a 2D output file, it makes a
# 3D surface plot.
#
# Usage:
#   plot_result_3D.py file.out
# ------------------------------------------------------------


# first things first: check whether you can import the hydro python modules
from check_module_is_in_pythonpath import try_to_import

try_to_import()


from mesh_hydro_utils import get_only_cmdlinearg, label_to_kwargs
from mesh_hydro_io import read_output
from mesh_hydro_plotting import plot_1D, plot_2D_in_3D, save_plot

if __name__ == "__main__":

    fname = get_only_cmdlinearg()
    ndim, rho, u, p, t, step = read_output(fname)

    if ndim == 1:
        fig = plot_1D(rho, u, p)
        save_plot(fig, fname)
    elif ndim == 2:
        fig = plot_2D_in_3D(rho, u, p, t=t)
        save_plot(fig, fname, case="3D")
