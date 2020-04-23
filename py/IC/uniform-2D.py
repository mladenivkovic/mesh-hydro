#!/usr/bin/env python3

#---------------------------------------------------
# Create 2D uniform ICs. Everything is identical.
#---------------------------------------------------


import numpy as np
from hydro_io import write_ic


nx = 100
rho_all = 1.
ux_all = 0.
uy_all = 0.
p_all = 1.




rho = np.ones((nx, nx), dtype=np.float) * rho_all
u = np.zeros((nx, nx, 2), dtype=np.float)
u[:, :, 0] = ux_all
u[:, :, 1] = uy_all
p = np.ones((nx, nx), dtype=np.float) * p_all



write_ic("uniform-2D-{0:d}.dat".format(nx), 2, rho, u, p)
