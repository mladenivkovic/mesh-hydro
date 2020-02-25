#!/usr/bin/env python3



#----------------------------------------------------
# Plot results of a 1D or 2D file as 1D solutions,
# and overplot the analytical solution to the Riemann
# problem.
#
# Usage:
#   ./plot_riemann_result.py <outputfile> <IC-file>
#
# it needs the IC file as well in order do be able to
# solve the riemann problem. If you just want to plot
# the results, use the appropriate scripts.
#----------------------------------------------------



from hydro_io import read_ic, read_output, check_file_exists
from hydro_riemann import riemann_solver
from hydro_plotting import plot_1D
from sys import argv



if __name__ == "__main__":
    
    fname = argv[1]
    check_file_exists(fname)
    icfname = argv[2]
    check_file_exists(icfname)


    ndim, rho, u, p, t, step = read_output(fname)

    nsim, twostate, rhoIC, uIC, pIC = read_ic(icfname, nx=rho.shape[0])

    if twostate:
        fig = plot_1D(rho, u, p, fname, t=t, nosave=True, draw_legend=True)

        rho_sol, u_sol, p_sol = riemann_solver(rhoIC, uIC, pIC, t)
        kwargs = { "linestyle":"--", }
        plot_1D(rho_sol, u_sol, p_sol, fname, t="python solver", nosave=False, draw_legend=True, fig=fig, kwargs = kwargs)

    else:
        print("Can't work with non-riemann ICs.")
