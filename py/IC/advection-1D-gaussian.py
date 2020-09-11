#!/usr/bin/env python3

# ---------------------------------------------------
# Create IC conditions for 1D advection where the
# profile is a gaussian
# ---------------------------------------------------


import numpy as np
from mesh_hydro_io import write_ic


nx = 100


rho = np.zeros(nx, dtype=np.float)
u = np.ones(nx, dtype=np.float)
p = np.ones(nx, dtype=np.float)


dx = 1.0 / nx

for i in range(nx):
    x = (i + 0.5) * dx
    rho[i] = 1 + 2 * np.exp(-((x - 0.5) ** 2 / 0.05))


write_ic("advection-1D-gaussian-{0:d}.dat".format(nx), 1, rho, u, p)
