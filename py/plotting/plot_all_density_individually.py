#!/usr/bin/env python3


#------------------------------------------------------------------------------------
# Create a plot of the result for every file given as cmdline argument.
# Plots all results in individual images.
#
# Usage:
#   plot_all_density_individually.py <file1> <file2> ... <file N>
#------------------------------------------------------------------------------------



from hydro_utils import get_only_cmdlinearg, get_all_files_with_same_basename
from hydro_io import read_output
from hydro_plotting import plot_1D_density_only, plot_2D_density_only

from sys import argv


# plotting parameters
dots = False    # overplot dots on 1D plot


if __name__ == "__main__":
    
    filelist = argv[1:]

    for fname in filelist:

        ndim, rho, u, p, t = read_output(fname)

        if ndim == 1:
            plot_1D_density_only(rho, fname, dots=dots)
        elif ndim == 2:
            plot_2D_density_only(rho, fname)
