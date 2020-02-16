#!/usr/bin/env python3


#-------------------------------------
# Create a plot of an output file.
# Plots only the density.
# Usage:
#   plot_density.py file.out
#-------------------------------------



from hydro_utils import get_only_cmdlinearg
from hydro_io import read_output
from hydro_plotting import plot_1D_density_only, plot_2D_density_only


# plotting parameters
dots = False    # overplot dots on 1D plot


if __name__ == "__main__":
    
    fname = get_only_cmdlinearg()
    ndim, rho, u, p, t = read_output(fname)

    if ndim == 1:
        plot_1D_density_only(rho, fname, dots=dots, t=t, draw_legend=True)
    elif ndim == 2:
        plot_2D_density_only(rho, fname, t=t)
