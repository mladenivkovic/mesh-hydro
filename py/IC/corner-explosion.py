#!/usr/bin/env python3

# ---------------------------------------------------
# Create IC conditions for an explosion in the lower
# left corner of the box.
# ---------------------------------------------------


import numpy as np
from mesh_hydro_io import write_ic


nx = 200


rho = np.ones((nx, nx), dtype=np.float)
u = np.zeros((nx, nx, 2), dtype=np.float)
p = np.ones((nx, nx), dtype=np.float) * 1e-5

dx = 1.0 / nx

p[0, 0] = 1.0 / dx

write_ic("corner-explosion-{0:d}.dat".format(nx), 2, rho, u, p)
