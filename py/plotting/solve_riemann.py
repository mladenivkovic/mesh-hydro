#!/usr/bin/env python3



#--------------------------------------------------------
# Solve Riemann problem and plot results of a 1D 
# two-state IC file and sample the solution at given time t
#
# Usage:
#   ./solve_riemann.py <ic file> <t>
#--------------------------------------------------------



from hydro_io import read_ic, read_output, check_file_exists
from hydro_riemann import riemann_solver
from hydro_plotting import plot_1D
from sys import argv
from matplotlib import pyplot as plt



if __name__ == "__main__":
    
    fname = argv[1]
    check_file_exists(fname)
    t = float(argv[2])


    nsim, twostate, rhoIC, uIC, pIC = read_ic(fname, nx=200)

    if twostate:
        rho_sol, u_sol, p_sol = riemann_solver(rhoIC, uIC, pIC, t)
        fig = plot_1D(rho_sol, u_sol, p_sol, fname, t=t, nosave=True, draw_legend=True,)
        plt.savefig(fname[:-3]+"png")

    else:
        print("Can't work with non-riemann ICs.")
