#!/usr/bin/env python3

#---------------------------------------------------
# Create IC conditions for an explosion in the lower
# left corner of the box.
#---------------------------------------------------


import numpy as np
from hydro_io import write_ic


nx = 200


rho = np.ones((nx, nx), dtype=np.float)
u = np.zeros((nx, nx, 2), dtype=np.float)
p = np.ones((nx, nx), dtype=np.float)*1e-5

dx = 1./nx

if nx % 2 == 0:
    c = int(nx / 2) - 1
    p[c,c] = 1./dx
    p[c+1, c] = 1./dx
    p[c,c+1] = 1./dx
    p[c+1,c+1] = 1./dx
else:
    c = int(nx / 2)
    p[c, c] = 1./dx
print("c=", c, "nx", nx)

write_ic("center-explosion-{0:d}.dat".format(nx), 2, rho, u, p)
