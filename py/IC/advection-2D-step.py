#!/usr/bin/env python3

# ---------------------------------------------------
# Create IC conditions for 2D advection where the
# profile is just a step function
# ---------------------------------------------------


import numpy as np
from mesh_hydro_io import write_ic


nx = 100


rho = np.ones((nx, nx), dtype=np.float)
u = np.zeros((nx, nx, 2), dtype=np.float)
p = np.ones((nx, nx), dtype=np.float)

u[:, :, 0] = 0.2
u[:, :, 1] = 0.2


dx = 1.0 / nx

for i in range(nx):
    xcenter = (i + 0.5) * dx
    if xcenter < 1.0 / 3:
        continue
    elif xcenter > 2.0 / 3:
        break
    for j in range(nx):
        ycenter = (j + 0.5) * dx
        if ycenter < 1.0 / 3:
            continue
        elif ycenter < 2.0 / 3:
            rho[i, j] = 2
        else:
            break

write_ic("advection-2D-step-{0:d}.dat".format(nx), 2, rho, u, p)
