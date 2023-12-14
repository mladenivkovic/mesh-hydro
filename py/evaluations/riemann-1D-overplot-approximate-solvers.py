#!/usr/bin/env python3


# ------------------------------------------------------------------------------------
# Plot multiple approximate solvers and exact Riemann solver solution into one
# plot to compare.
# Needs the initial conditions file, and builds the filenames for the differnt
# solvers based on the IC file name.
# ------------------------------------------------------------------------------------

# first things first: check whether you can import the hydro python modules
from check_module_is_in_pythonpath import try_to_import

try_to_import()


from mesh_hydro_io import read_output
from mesh_hydro_plotting import plot_1D, save_plot
from mesh_hydro_utils import label_to_kwargs

from sys import argv


icfile = argv[1]


solvers = ["EXACT", "TRRS", "TSRS"]
file_prefix = icfile[:-4]
filelist = [file_prefix + "-" + s + "-0001.out" for s in solvers]
labels = solvers
linestyles = ["-", "--", ":"]


if __name__ == "__main__":
    fig = None
    for i, f in enumerate(filelist):
        ndim, rho, u, p, t, step = read_output(f)

        if ndim != 1:
            print("I can't overplot 2D stuff...")
            quit(1)
        else:
            kwargs = label_to_kwargs(labels[i])
            kwargs["linestyle"] = linestyles[i]
            fig = plot_1D(rho, u, p, draw_legend=True, fig=fig, kwargs=kwargs)

    fig.suptitle("t = {0:.3f}".format(t))
    figname = "riemann-approximate-" + file_prefix + ".png"
    save_plot(fig, f, fname_force=figname)
