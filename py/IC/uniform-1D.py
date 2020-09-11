#!/usr/bin/env python3

# ---------------------------------------------------
# Create 1D uniform ICs. Everything is identical.
# ---------------------------------------------------


import numpy as np
from mesh_hydro_io import write_ic


nx = 100
rho_all = 1.0
ux_all = 0.0
p_all = 1.0


rho = np.ones(nx, dtype=np.float) * rho_all
u = np.ones(nx, dtype=np.float) * ux_all
p = np.ones(nx, dtype=np.float) * p_all


write_ic("uniform-1D-{0:d}.dat".format(nx), 1, rho, u, p)
