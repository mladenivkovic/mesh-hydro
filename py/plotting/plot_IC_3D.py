#!/usr/bin/env python3


#------------------------------------------------------------
# Create a plot of an output file.
# If it's a 1D output file, it creates the same as the 
# plot_result.py script. If it's a 2D output file, it makes a
# 3D surface plot.
# 
# Usage:
#   plot_result_3D.py file.out
#-------------------------------------



from hydro_utils import get_only_cmdlinearg
from hydro_io import read_ic
from hydro_plotting import plot_1D, plot_2D_in_3D


# plotting parameters
dots = False    # overplot dots on 1D plot


if __name__ == "__main__":
    
    fname = get_only_cmdlinearg()
    ndim, twostate, rho, u, p = read_ic(fname)

    if twostate:
        plot_1D(rho, u, p, fname)
    else:
        if ndim == 1:
            plot_1D(rho, u, p, fname, dots=dots)
        elif ndim == 2:
            plot_2D_in_3D(rho, u, p, fname)
