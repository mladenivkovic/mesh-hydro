#!/usr/bin/env python3


# -----------------------------------------------------------------
# Get a plot of a Riemann problem with all 3 waves present,
# and mark the waves appropriately
# -----------------------------------------------------------------


from hydro_riemann import riemann_solver, soundspeed, find_star_state
from hydro_plotting import plot_1D, save_plot
import numpy as np
from matplotlib import pyplot as plt


# --------------------------
# setup initial conditions
# --------------------------

rhoL = 1.0
uL = 0.0
pL = 1.0

rhoR = 0.125
uR = 0.0
pR = 0.1

t = 0.20

nx = 500


rho = np.empty((nx), dtype=np.float)
u = np.empty((nx), dtype=np.float)
p = np.empty((nx), dtype=np.float)

nxhalf = nx // 2
rho[:nxhalf] = rhoL
u[:nxhalf] = uL
p[:nxhalf] = pL
rho[nxhalf:] = rhoR
u[nxhalf:] = uR
p[nxhalf:] = pR


# solve the Riemann problem

rho_sol, u_sol, p_sol = riemann_solver(rho, u, p, t)

# find minima and maxima values for plots
mins = []
maxs = []

for arr in [rho_sol, u_sol, p_sol]:
    mins.append(min(0.95 * arr.min(), -0.05 * arr.max()))
    maxs.append(1.05 * arr.max())


# --------------------------
# Now for the plotting
# --------------------------


# First plot the Initial Conditions
fig = plot_1D(rho, u, p)

for i, ax in enumerate(fig.axes):
    ax.set_xlim(0, 1)
    ax.set_ylim(mins[i], maxs[i])
fig.suptitle("Initial Conditions")

save_plot(fig, fname_force="riemann_IC.pdf")
plt.close()


# Now plot the Solution
fig = plot_1D(rho_sol, u_sol, p_sol)
fig.suptitle("Solution of Riemann problem at t = {0:.2f}".format(t))


# add lines for waves

pstar, ustar = find_star_state(rhoL, uL, pL, rhoR, uR, pR)
gamma = 5.0 / 3.0
alpha = 0.5 * (gamma - 1) / gamma

aL = soundspeed(pL, rhoL)
aR = soundspeed(pR, rhoR)
astarL = aL * (pstar / pL) ** alpha

SHL = uL - aL  # speed of the head of left rarefaction
STL = ustar - astarL  # speed of the tail of left rarefaction
SR = uR + aR * np.sqrt(
    0.5 * (gamma + 1) / gamma * pstar / pR + alpha
)  # right shock speed


# find positions of waves as given time. Add +0.5 because they originate at x = 0.5
xrarehead = 0.5 + SHL * t
xraretail = 0.5 + STL * t
xcontact = 0.5 + ustar * t
xshock = 0.5 + SR * t


for i, ax in enumerate(fig.axes):
    ax.set_xlim(0, 1)
    ax.set_ylim(mins[i], maxs[i])
    ax.plot([xrarehead, xrarehead], [mins[i], maxs[i]], "r:", lw=1)
    ax.plot([xraretail, xraretail], [mins[i], maxs[i]], "r:", lw=1)
    ax.plot([xcontact, xcontact], [mins[i], maxs[i]], "g:", lw=1)
    ax.plot([xshock, xshock], [mins[i], maxs[i]], ":", color="orange", lw=1)


save_plot(fig, fname_force="riemann_exact_solution.pdf")
plt.close()
