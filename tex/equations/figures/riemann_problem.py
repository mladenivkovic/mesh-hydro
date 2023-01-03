#!/usr/bin/env python3

# -----------------------------------------------------------------
# Create a plot for the Riemann IVP
# -----------------------------------------------------------------


import numpy as np
from matplotlib import pyplot as plt
import matplotlib

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
    "figure.subplot.left": 0.05,
    "figure.subplot.right": 0.97,
    "figure.subplot.bottom": 0.12,
    "figure.subplot.top": 0.92,
    "figure.subplot.wspace": 0.25,
    "figure.subplot.hspace": 0.25,
    "figure.dpi": 200,
    "lines.markersize": 6,
    "lines.linewidth": 2.0,
}

matplotlib.rcParams.update(params)


fig = plt.figure(figsize=(6, 3))
ax = fig.add_subplot(111)
ax.plot([-10, 10], [0, 0], "k")
ax.plot([0, 0], [0, 10], "k", lw=1)

ax.axis("off")
plt.text(0, -1, r"0", verticalalignment="center", horizontalalignment="center")
plt.text(10.5, 0.0, r"x", verticalalignment="center", horizontalalignment="center")
plt.text(
    -5, 5, r"$\mathbf{U}_L$", horizontalalignment="center", verticalalignment="center"
)
plt.text(
    5, 5, r"$\mathbf{U}_R$", horizontalalignment="center", verticalalignment="center"
)


plt.savefig("riemann_problem.pdf")
