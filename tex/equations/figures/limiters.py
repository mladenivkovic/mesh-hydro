#!/usr/bin/env python3

#----------------------------------------------
# Create r, phi(r) plots for various limiters
#----------------------------------------------



import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.patches as patches

params = {
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'font.size': 12,
    'font.family': 'serif',
    'font.serif': 'DejaVu Sans',
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'figure.subplot.left'    : 0.10,
    'figure.subplot.right'   : 0.97,
    'figure.subplot.bottom'  : 0.12,
    'figure.subplot.top'     : 0.92,
    'figure.subplot.wspace'  : 0.25,
    'figure.subplot.hspace'  : 0.25,
    'figure.dpi' : 200,
    'lines.markersize' : 6,
    'lines.linewidth' : 2.
}

matplotlib.rcParams.update(params)





def minmod(r):
    # compute minmod(1, r)

    if 1*r <= 0:
        return 0
    if 1 < r:
        return 1
    if r < 1:
        return r


def superbee(r):
    # get superbee limiter
    return max(0, min(1, 2*r), min(2, r))


def MC(r):
    # get monotinized central-difference
    return max(0, min(0.5*(1+r), 2, 2*r))


def vanleer(r):
    # get van Leer limiter
    return (r + abs(r))/(1 + abs(r))




r = np.linspace(-1, 5, 1000)

fig = plt.figure()
ax = fig.add_subplot(111)


# draw allowable regions

ax.set_xlim(-1, 5)
ax.set_ylim(-0.2, 2.2)
ax.add_patch(patches.Polygon([(0, 0), (1, 1), (0.5, 1)], fill=True, facecolor='darkgrey'))
ax.add_patch(patches.Polygon([(1, 1), (6, 1), (6, 2), (2, 2)], fill=True, facecolor='darkgrey'))

# draw asymptotes
ax.plot([0, 1.1], [0, 2.2], 'k--', lw=3)
ax.annotate(r"$\phi(r)=2r$", [-0.05, 1.5])
ax.plot([0, 2.2], [0, 2.2], 'k--', lw=3)
ax.annotate(r"$\phi(r)=r$", [0.7, 0.5])
ax.plot([-1, 5], [1, 1], 'k--', lw=3)
ax.annotate(r"$\phi(r)=1$", [-0.9, 0.9])
ax.plot([-1, 5], [2, 2], 'k--', lw=3)
ax.annotate(r"$\phi(r)=2$", [-0.9, 2.05])



for name, phi in zip(["minmod", "superbee", "mc", "van Leer"], [minmod, superbee, MC, vanleer]):
    y = [phi(x) for x in r]
    ax.plot(r, y, label=name)

ax.legend(loc="lower right")

ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$\phi(r)$")

plt.savefig("limiters.pdf")
