#!/usr/bin/env python3

# -----------------------------------------------------------------
# Create a plot for the solution structure of the Riemann IVP
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

xmin = -10
xmax = -xmin
ymin = 0
ymax = 10

fig = plt.figure(figsize=(6, 3))
ax = fig.add_subplot(111)
ax.plot([xmin, xmax], [0, 0], "k", zorder=10)
ax.plot([0, 0], [ymin, ymax], "k", zorder=10)


def x_for_ymax(a):
    res = ymax / a
    if res > xmax:
        return xmax
    elif res < xmin:
        return xmin
    else:
        return res


# contact wave
ax.plot([0, x_for_ymax(4)], [0, ymax], "b")
# right wave
ax.plot([0, x_for_ymax(0.9)], [0, ymax], "r")
# left wave
ax.plot([0, x_for_ymax(-0.9)], [0, ymax], "r")

ax.axis("off")
# axis
# origin
plt.text(0, -1, r"0", verticalalignment="center", horizontalalignment="center")
# x label
plt.text(10.5, 0.0, r"x", verticalalignment="center", horizontalalignment="center")
# t label
plt.text(0, 11, r"t", verticalalignment="center", horizontalalignment="center")

# U_L
plt.text(
    -6, 3, r"$\mathbf{U}_L$", horizontalalignment="center", verticalalignment="center"
)
# U_R
plt.text(
    6, 3, r"$\mathbf{U}_R$", horizontalalignment="center", verticalalignment="center"
)
# U*_L
plt.text(
    -2, 6, r"$\mathbf{U}_L^*$", horizontalalignment="center", verticalalignment="center"
)
# U*_R
plt.text(
    4, 6, r"$\mathbf{U}_R^*$", horizontalalignment="center", verticalalignment="center"
)

# wave 1
plt.text(
    -9,
    7.5,
    r"$(1)$",
    color="r",
    horizontalalignment="center",
    verticalalignment="center",
)
# wave 2
plt.text(
    1.5,
    9,
    r"$(2)$",
    color="b",
    horizontalalignment="center",
    verticalalignment="center",
)
# wave 3
plt.text(
    9,
    7.5,
    r"$(3)$",
    color="r",
    horizontalalignment="center",
    verticalalignment="center",
)


plt.savefig("riemann_solution.pdf")
