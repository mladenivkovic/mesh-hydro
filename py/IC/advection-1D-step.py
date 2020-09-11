#!/usr/bin/env python3

# ---------------------------------------------------
# Create IC conditions for 1D advection where the
# profile is just a step function
# ---------------------------------------------------


import numpy as np
from mesh_hydro_io import write_ic


nx = 20000


rho = np.ones(nx, dtype=np.float)
u = np.ones(nx, dtype=np.float)
p = np.ones(nx, dtype=np.float)


dx = 1.0 / nx

for i in range(nx):
    center = (i + 0.5) * dx
    if center < 1.0 / 3:
        continue
    elif center < 2.0 / 3:
        rho[i] *= 2
    else:
        break


write_ic("advection-1D-step-{0:d}.dat".format(nx), 1, rho, u, p)
