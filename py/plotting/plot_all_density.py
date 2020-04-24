#!/usr/bin/env python3


#------------------------------------------------------------------------------------
# Create a plot of all output files with the same base name as given file in 
# this directory. Basename in this case means everything before _XXXX.out
# If more than one cmdline arg is given, it will instead interpret every arg as an
# individual file, and only plot those.
#
# If the first two files have the same output time, the lines will be labelled by
# the file name, not the output time.
#
# ONLY WORKS FOR 1D OUTPUT.
# Plots all results in one plot.
#
# Usage:
#   plot_all_results.py file.out     # will find files with same basename
# or:
#   plot_all_results.py <file1> <file2> ... <file N>
#------------------------------------------------------------------------------------


# first things first: check whether you can import the hydro python modules
from check_module_is_in_pythonpath import try_to_import
try_to_import()


from hydro_utils import get_only_cmdlinearg, get_all_files_with_same_basename, label_to_kwargs
from hydro_io import read_output
from hydro_plotting import plot_1D_density_only, save_plot

from sys import argv


if __name__ == "__main__":
    
    if len(argv) == 2:
        fname = get_only_cmdlinearg()
        filelist = get_all_files_with_same_basename(fname)
    else:
        filelist = argv[1:]


    fig = None
    tlist = []
    label_is_fname = False
    i = 0
    while i < len(filelist):
        f = filelist[i]

        ndim, rho, u, p, t, step = read_output(f)

        # store first time and file name
        if not label_is_fname:
            if t not in tlist:
                tlist.append(t)
                labelval = t
            else:
                label_is_fname = True
                labelval = f
                fig = None
                i = 0 # restart!
                continue
        else:
            labelval = f

        if ndim != 1:
            print("I can't overplot 2D stuff...")
            quit(1)
        else:
            fig = plot_1D_density_only(rho, draw_legend=True, fig = fig, kwargs=label_to_kwargs(labelval), )

    save_plot(fig, f, case="density-overplotted")
