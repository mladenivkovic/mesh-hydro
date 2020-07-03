#!/usr/bin/env python3



#--------------------------------------------------------------------------------
# Plot results of a 1D or 2D file as 1D solutions, and overplot the analytical 
# solution to the Riemann problem for an arbitrary number of files.
#
# Usage:
#   ./plot_all_riemann_result.py <file1> <file2> ... <fileN> <IC-file>
#
# it needs the IC file as well in order do be able to solve the riemann problem. 
# If you just want to plot the results, use the appropriate scripts.
#--------------------------------------------------------------------------------


# first things first: check whether you can import the hydro python modules
from check_module_is_in_pythonpath import try_to_import
try_to_import()



from mesh_hydro_io import read_ic, read_output, check_file_exists
from mesh_hydro_utils import label_to_kwargs
from mesh_hydro_riemann import riemann_solver
from mesh_hydro_plotting import plot_1D, save_plot
from sys import argv
import numpy as np



if __name__ == "__main__":
    
    # last file needs to be IC file
    nfiles = len(argv) - 2
    icfname = argv[-1]

    fig = None

    for fname in argv[1:-1]:
        check_file_exists(fname)
        ndim, rho, u, p, t, step = read_output(fname)

        label = fname.replace("_", "-")

        if ndim == 2:
            rho = np.mean(rho, axis=0)
            u = np.mean(u[:, :, 0], axis=0)
            p = np.mean(p, axis = 0)
            label = fname.replace("_", "-") + r"; mean value along y"

        kwargs = label_to_kwargs(label)
        fig = plot_1D(rho, u, p, draw_legend=True, kwargs=kwargs, fig=fig)



    check_file_exists(icfname)
    ndim, twostate, rhoIC, uIC, pIC = read_ic(icfname, nx=rho.shape[0])

    if twostate:
        rho_sol, u_sol, p_sol = riemann_solver(rhoIC, uIC, pIC, t)
        kwargs = { "linestyle":"--", }
        kwargs = label_to_kwargs(t="python solver", kwargs=kwargs)
        fig = plot_1D(rho_sol, u_sol, p_sol, draw_legend=True, fig=fig, kwargs = kwargs)

        save_plot(fig, fname, case='overplotted')

    else:
        print("Can't work with non-riemann ICs.")
