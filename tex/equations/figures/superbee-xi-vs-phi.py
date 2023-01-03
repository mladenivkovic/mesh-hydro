#!/usr/bin/env python3

from matplotlib import pyplot as plt
import matplotlib
import numpy as np

# -------------------------------------------
# Plot phi(r) and xi(r) w.r.t. r for
# different values of omega
# -------------------------------------------


params = {
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "font.size": 14,
    "font.family": "serif",
    "font.serif": "DejaVu Sans",
    "legend.fontsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "text.usetex": True,
    "figure.subplot.left": 0.15,
    "figure.subplot.right": 0.95,
    "figure.subplot.bottom": 0.12,
    "figure.subplot.top": 0.92,
    "figure.subplot.wspace": 0.25,
    "figure.subplot.hspace": 0.25,
    "figure.dpi": 200,
    "lines.markersize": 6,
    "lines.linewidth": 2.0,
}
matplotlib.rcParams.update(params)


def superbee_xi(r, omega):

    xi = np.zeros(r.shape)
    xi[r > 0] = 2 * r[r > 0]
    xi[r > 0.5] = 1.0

    denor = 1.0 - omega + (1.0 + omega) * r
    xiR = 2.0 / denor
    val = np.minimum(xiR, r)
    val = np.minimum(val, 2.0)
    xi[r > 1] = val[r > 1]

    return xi


def superbee_phi(r):

    n = r.shape[0]

    phi = np.maximum(0, np.minimum(np.ones(n), 2.0 * r))
    phi = np.maximum(phi, np.minimum(2.0 * np.ones(n), r))
    return phi


r = np.linspace(-1, 5, 1000)

phi = superbee_phi(r)
xi_upwind = superbee_xi(r, 1)
xi_centered = superbee_xi(r, 0)
xi_downwind = superbee_xi(r, -1)


plt.figure()
plt.plot(r, phi, label=r"$\phi(r)$")
plt.plot(r, xi_upwind, label=r"$\xi(r)$, $\omega = 1$")
plt.plot(r, xi_centered, label=r"$\xi(r)$, $\omega = 0$")
plt.plot(r, xi_downwind, label=r"$\xi(r)$, $\omega = -1$")
plt.legend()

plt.xlabel(r"$r$")
plt.ylabel("[1]")
#  plt.show()
plt.savefig("superbee-xi-vs-phi.pdf")
