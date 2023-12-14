#!/usr/bin/env python3


# --------------------------------------------------------
# Solve Riemann problem and plot results of a 1D
# two-state IC file and sample the solution at given time t
#
# Usage:
#   ./solve_riemann.py <ic file> <t>
# --------------------------------------------------------


# first things first: check whether you can import the hydro python modules
from check_module_is_in_pythonpath import try_to_import

try_to_import()


from mesh_hydro_io import read_ic, read_output, check_file_exists
from mesh_hydro_riemann import riemann_solver
from mesh_hydro_plotting import plot_1D, save_plot
from mesh_hydro_utils import label_to_kwargs
from sys import argv
from matplotlib import pyplot as plt


if __name__ == "__main__":
    fname = argv[1]
    check_file_exists(fname)
    t = float(argv[2])

    nsim, twostate, rhoIC, uIC, pIC = read_ic(fname, nx=200)

    if twostate:
        rho_sol, u_sol, p_sol = riemann_solver(rhoIC, uIC, pIC, t)
        fig = plot_1D(
            rho_sol, u_sol, p_sol, draw_legend=True, kwargs=label_to_kwargs(t)
        )
        save_plot(fig, fname, case="riemann-solver")

    else:
        print("Can't work with non-riemann ICs.")
