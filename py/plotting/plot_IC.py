#!/usr/bin/env python3


# -------------------------------------
# Create a plot of an IC file.
# Usage:
#   plot_IC.py file.dat
# -------------------------------------


# first things first: check whether you can import the hydro python modules
from check_module_is_in_pythonpath import try_to_import

try_to_import()


from mesh_hydro_utils import get_only_cmdlinearg, label_to_kwargs
from mesh_hydro_io import read_ic
from mesh_hydro_plotting import plot_1D, plot_2D, save_plot


if __name__ == "__main__":

    fname = get_only_cmdlinearg()
    ndim, twostate, rho, u, p = read_ic(fname)

    if ndim == 1:
        fig = plot_1D(rho, u, p,)
    elif ndim == 2:
        fig = plot_2D(rho, u, p,)

    save_plot(fig, fname)
