#!/usr/bin/env python3


#----------------------------------------------------
# Create a plot containing different Riemann solver
# solutions with appropriate labels.
#
# expects 1 cmdline arg: the prefix of the IC file,
# from which the riemann program will generate a 
# file name.
#----------------------------------------------------



from hydro_io import read_ic, read_output, check_file_exists
from hydro_utils import label_to_kwargs
from hydro_riemann import riemann_solver
from hydro_plotting import plot_1D, save_plot
from sys import argv
import numpy as np
import os



def get_all_files_with_same_basename(fname):
    """
    Get a list of all files with the same basename as given file <fname>.
    Basename in this case means everything before -RIEMANN-ABCD-XXXX.out, which is the
    format of the riemann output files.
    also get a list of which solvers are used, which is stored in ABCD of the filename format
    described above.
    
    Returns: [filelist], [solvername list]
    """

    
    basedir = os.path.dirname(fname)
    if basedir == '':
        basedir = os.getcwd()

    allfiles = os.listdir(basedir)

    filelist = []

    for f in allfiles:
        if f.startswith(fname+"-RIEMANN") and f.endswith("0001.out"):
            filelist.append(f)

    filelist.sort()

    start = len(fname)+len("-RIEMANN-")
    end = len("-0001.out")
    solvernamelist = [f[start:-end] for f in filelist]

    return filelist, solvernamelist





if __name__ == "__main__":
    
    icprefix = argv[1]
    icfname = os.path.join("IC", icprefix+".dat")

    filelist, namelist = get_all_files_with_same_basename(icprefix)

    fig = None

    for i in range(len(filelist)):
        fname = filelist[i]
        solver = namelist[i]

        ndim, rho, u, p, t, step = read_output(fname)
        kwargs = label_to_kwargs(solver)
        fig = plot_1D(rho, u, p, draw_legend=True, fig=fig, kwargs=kwargs)




    # plot exact python solution
    nsim, twostate, rhoIC, uIC, pIC = read_ic(icfname, nx=rho.shape[0])

    if twostate:

        rho_sol, u_sol, p_sol = riemann_solver(rhoIC, uIC, pIC, t)
        kwargs = { "linestyle":"--", "color":'k'}
        kwargs = label_to_kwargs(t="python solver", kwargs=kwargs)
        fig = plot_1D(rho_sol, u_sol, p_sol, draw_legend=True, fig=fig, kwargs = kwargs)

        save_plot(fig, fname_force = icprefix+".png")

    else:
        print("Can't work with non-riemann ICs.")
