#!/usr/bin/env python3


#----------------------------------------------------
# Create a plot containing different Riemann solver
# solutions with appropriate labels.
#
# expects 2 cmdline args: the prefix of the IC file,
# from which the riemann program will generate a 
# file name; and the used solver which appears in
# the file name, e.g. GODUNOV
#----------------------------------------------------



from hydro_io import read_ic, read_output, check_file_exists
from hydro_utils import label_to_kwargs
from hydro_riemann import riemann_solver
from hydro_plotting import plot_1D, save_plot
from sys import argv
import numpy as np
import os



def get_all_files_with_same_basename(fname, solver):
    """
    Get a list of all files with the same basename as given file <fname>.
    Basename in this case means everything before <SOLVER>-ABCD-XXXX.out, which is the
    format of the riemann output files.
    also get a list of which RIEMANN solvers are used, (not hydro solvers! that should be the second arg)
    which is stored in ABCD of the filename format described above.
    
    Returns: [filelist], [solvername list]
    """

    
    basedir = os.path.dirname(fname)
    if basedir == '':
        basedir = os.getcwd()

    allfiles = os.listdir(basedir)

    filelist = []

    for f in allfiles:
        if f.startswith(fname+"-"+solver) and f.endswith("0001.out"):
            filelist.append(f)

    filelist.sort()

    start = len(fname)+len(solver) + 2 # +2 for 2 dashes between the solver and what comes after it
    end = len("-0001.out")
    riemannsolvernamelist = []
    for f in filelist:
        char = 'a'
        i = 0
        while char != '-':
            i += 1
            char = f[start + i]
        riemannsolvernamelist.append(f[start:start+i])
        # [:i] means i not included, but i shoud be a dash, so we're good


    return filelist, riemannsolvernamelist





if __name__ == "__main__":
    
    icprefix = argv[1]
    solver = argv[2]
    icfname = os.path.join("IC", icprefix+".dat")

    filelist, namelist = get_all_files_with_same_basename(icprefix, solver)

    fig = None

    linestyles = ["--", "-.", ":"]

    for i in range(len(filelist)):
        fname = filelist[i]
        solver = namelist[i]
        ls = linestyles[i % len(linestyles)]

        ndim, rho, u, p, t, step = read_output(fname)

        if ndim == 1:
            kwargs = label_to_kwargs(solver)
            kwargs["linestyle"] = ls
            fig = plot_1D(rho, u, p, draw_legend=True, fig=fig, kwargs=kwargs)
        elif ndim == 2:

            rho = np.mean(rho, axis=0)
            u = np.mean(u[:, :, 0], axis=0)
            p = np.mean(p, axis = 0)
            label = "t = {0:.3f}; mean value along y".format(t)
            kwargs = {"label":label}
            kwargs["linestyle"] = ls
            fig = plot_1D(rho, u, p, draw_legend=True, fig=fig, kwargs=kwargs)


    fig.suptitle(r"t = {0:.3f}".format(t))



    # plot exact python solution
    nsim, twostate, rhoIC, uIC, pIC = read_ic(icfname, nx=rho.shape[0])

    if twostate:

        rho_sol, u_sol, p_sol = riemann_solver(rhoIC, uIC, pIC, t)
        kwargs = { "linestyle":"-", "color":'k', "zorder":-1}
        kwargs = label_to_kwargs(t="python solver", kwargs=kwargs)
        fig = plot_1D(rho_sol, u_sol, p_sol, draw_legend=True, fig=fig, kwargs = kwargs)

        save_plot(fig, fname_force = "GODUNOV-"+icprefix+"-{0:1d}D.png".format(ndim))

    else:
        print("Can't work with non-riemann ICs.")
