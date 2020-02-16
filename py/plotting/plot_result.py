#!/usr/bin/env python3


#-------------------------------------
# Create a plot of an output file.
# Usage:
#   plot_result.py file.out
#-------------------------------------



from hydro_utils import get_only_cmdlinearg
from hydro_io import read_output
from hydro_plotting import plot_1D, plot_2D


# plotting parameters
dots = False    # overplot dots on 1D plot


if __name__ == "__main__":
    
    fname = get_only_cmdlinearg()
    ndim, rho, u, p, t = read_output(fname)


    if ndim == 1:
        plot_1D(rho, u, p, fname, dots=dots)
    elif ndim == 2:
        plot_2D(rho, u, p, fname)
