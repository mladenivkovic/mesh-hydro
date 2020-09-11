#!/usr/bin/env python3

# ---------------------------------------------------
# Create IC conditions for 1D advection where the
# profile contains a step function, a triangle, a
# gaussian and a sine
# ---------------------------------------------------


import numpy as np
from mesh_hydro_io import write_ic


nx = 100

baseline = 1  # baseline density
amplitude = 1.5  # initial density peak


rho = np.ones(nx, dtype=np.float)
u = np.ones(nx, dtype=np.float)
p = np.ones(nx, dtype=np.float)


dx = 1.0 / nx

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

for i in range(nx):

    if i - region * width >= width:
        region += 1

    if region == 0:
        # empty
        continue

    elif region == 1:
        # step function
        rho[i] = amplitude

    elif region == 2:
        # empty
        continue

    elif region == 3:
        # triangle
        center = (region + 0.5) * width  # center index
        if i <= center:
            a = (amplitude - baseline) / (width * dx / 2)
            b = baseline - a * (region * width + 0.5) * dx
            x = (i + 0.5) * dx
            rho[i] = a * x + b

        else:
            a = (baseline - amplitude) / (width * dx / 2)
            b = baseline - a * ((region + 1) * width + 0.5) * dx
            x = (i + 0.5) * dx
            rho[i] = a * x + b

    elif region == 4:
        # empty
        continue

    elif region == 5:
        # gauss
        center = (
            (region + 0.5) * width + 0.5
        ) * dx  # value of x at the center of this shape
        x = (i + 0.5) * dx
        rho[i] = baseline + (amplitude - baseline) * np.exp(
            -((x - center) ** 2 / 0.0007)
        )

    elif region == 6:
        # empty
        continue

    elif region == 7:
        # sin region.
        rho[i] = baseline + (amplitude - baseline) * np.sin(
            (i - region * width + 0.5) / width * np.pi
        )
        continue

    else:
        # empty
        continue


write_ic("advection-1D-four-shapes-{0:d}.dat".format(nx), 1, rho, u, p)
