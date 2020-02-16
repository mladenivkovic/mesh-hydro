#!/usr/bin/env python3

#---------------------------------------------------
# Create IC conditions for 1D advection where the
# profile contains a step function, a triangle, a
# gaussian and a sine
#---------------------------------------------------


import numpy as np
from hydro_io import write_ic


nx = 200

baseline = 1.   # baseline density
amplitude = 1.5 # initial density peak



rho = np.ones((nx, nx), dtype=np.float)
u = np.ones((nx, nx, 2), dtype=np.float)
p = np.ones((nx, nx), dtype=np.float)



dx = 1./nx

# we add 8 regions of equal width:
# 1) empty
# 2) step function
# 3) empty
# 4) triangle
# 5) gauss
# 6) empty
# 7) sin
region = 0

width = nx // 8

region_i = 0
for i in range(nx):

    if i - region_i * width >= width:
        region_i += 1

    region_j = 0

    for j in range(nx):

        if j - region_j * width >= width:
            region_j += 1

        if region_i != region_j:
            continue
        else:
            region = region_i

        if region == 0:
            # empty
            continue

        elif region == 1:
            # step function
            rho[i, j] = amplitude

        elif region == 2:
            # empty
            continue

        elif region == 3:
            # triangle
            center = ((region + 0.5) * width) # center index
            if i <= center:
                a = (amplitude - baseline) / (width*dx/2)
                b = baseline - a * (region * width + 0.5)*dx
                x = (i + 0.5)*dx
                rho[i, j] += (a*x + b - baseline) / 2

            else:
                a = (baseline - amplitude) / (width*dx/2)
                b = baseline - a * ((region+1) * width + 0.5)*dx
                x = (i + 0.5)*dx
                rho[i, j] += (a*x + b - baseline) / 2

            if j <= center:
                a = (amplitude - baseline) / (width*dx/2)
                b = baseline - a * (region * width + 0.5)*dx
                y = (j + 0.5)*dx
                rho[i, j] += (a*y + b - baseline) / 2

            else:
                a = (baseline - amplitude) / (width*dx/2)
                b = baseline - a * ((region+1) * width + 0.5)*dx
                y = (j + 0.5)*dx
                rho[i, j] += (a*y + b - baseline) / 2

        elif region == 4:
            # empty
            continue

        elif region == 5:
            # gauss
            center = ((region + 0.5) * width + 0.5) * dx # value of x at the center of this shape
            x = (i + 0.5) * dx
            y = (j + 0.5) * dx
            rho[i, j] = baseline + (amplitude - baseline) * np.exp(-((x-center)**2/0.0007)) * np.exp(-((y-center)**2/0.0007))

        elif region == 6:
            # empty
            continue

        elif region == 7:
            # sin region.
            rho[i, j] = baseline + (amplitude - baseline) * np.sin((i - region*width+0.5)/width*np.pi)* np.sin((j - region*width+0.5)/width*np.pi)
            continue

        else:
            # empty
            continue



write_ic("advection-2D-four-shapes.dat", 2, rho, u, p)
