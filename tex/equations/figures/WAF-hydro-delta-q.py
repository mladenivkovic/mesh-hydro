#!/usr/bin/env python3

# Create a plot for the delta q for the WAF method flux limiting


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
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'figure.subplot.left'    : 0.05,
    'figure.subplot.right'   : 0.95,
    'figure.subplot.bottom'  : 0.12,
    'figure.subplot.top'     : 0.92,
    'figure.subplot.wspace'  : 0.25,
    'figure.subplot.hspace'  : 0.25,
    'figure.dpi' : 200,
    'lines.markersize' : 6,
    'lines.linewidth' : 2.
}

matplotlib.rcParams.update(params)

xmin = -20
xmax = -xmin
dx = (xmax - xmin) * 0.01
ymin = 0
ymax = 10

fig = plt.figure(figsize=(12, 4))
ax = fig.add_subplot(111)
ax.axis("off")

ax.plot([xmin-dx, xmax+dx], [0, 0], 'k', zorder=10) # x axis
ax.plot([0, 0], [ymin, ymax+1], 'k', zorder=10)       # t axis



def x_for_ymax(a):
    # get correct x for a straight line given the slope a
    # cut off at [xmin, xmax]
    res = ymax/a
    if res > xmax:
        return xmax
    elif res < xmin:
        return xmin
    else:
        return res


def x_for_y(y, a):
    # get x for given y and slope a
    res = y/a
    return res


for origin, zaehler, sign in zip([0, -(xmax-xmin)/4, (xmax-xmin)/4], [1, 1, 3], ["+", "-", "+"]):

    # characteristics
    ax.plot([origin, origin+x_for_ymax(-4)], [0, ymax], 'k', lw=1)
    ax.plot([origin, origin+x_for_ymax(6)], [0, ymax], 'k', lw=1)
    ax.plot([origin, origin+x_for_ymax(2.5)], [0, ymax], 'k', lw=1)


    # dq1
    plt.text(origin-2, 11, 
        r"$\Delta q_{i"+sign+"{0:d}".format(zaehler)+r"/2}^{(1)}$",
        color='k',
        verticalalignment="top",
        horizontalalignment='center',)

    # dq2
    plt.text(origin+1.7, 11,
        r"$\Delta q_{i"+sign+"{0:d}".format(zaehler)+r"/2}^{(2)}$",
        color='k',
        verticalalignment="top",
        horizontalalignment='center',)

    # dq3
    plt.text(origin+4, 11,
        r"$\Delta q_{i"+sign+"{0:d}".format(zaehler)+r"/2}^{(3)}$",
        color='k',
        verticalalignment="top",
        horizontalalignment='center',)





# x label
plt.text(21, 0.0, r"$x$",
    verticalalignment="center",
    horizontalalignment='left',)
# t label
plt.text(0, 12, r"t",
    verticalalignment="center",
    horizontalalignment='center',)

dx = (xmax-xmin)/4
height = ymax - ymin

# U_i-1 label
plt.text(-1.5*dx, -1, r"$\mathbf{U}_{i-1}$",
    verticalalignment="center",
    horizontalalignment='center',)

# U_i label
plt.text(-.5*dx, -1, r"$\mathbf{U}_{i}$",
    verticalalignment="center",
    horizontalalignment='center',)

# U_i+1 label
plt.text(.5*dx, -1, r"$\mathbf{U}_{i+1}$",
    verticalalignment="center",
    horizontalalignment='center',)

# U_i+2 label
plt.text(1.5*dx, -1, r"$\mathbf{U}_{i+2}$",
    verticalalignment="center",
    horizontalalignment='center',)


# U_i-1 coloring
ax.add_patch(patches.Rectangle((xmin, ymin), width=dx, height=height, fill=True, facecolor='thistle'))
# U_i coloring
ax.add_patch(patches.Rectangle((xmin/2, ymin), width=dx, height=height, fill=True, facecolor='lightblue'))
# U_i+1 coloring
ax.add_patch(patches.Rectangle((0, ymin), width=(xmax-xmin)/4, height=(ymax-ymin), fill=True, facecolor='thistle'))
# U_i+2 coloring
ax.add_patch(patches.Rectangle((dx, ymin), width=dx, height=height, fill=True, facecolor='lightblue'))


#  plt.show()
plt.savefig("WAF-hydro-delta-q.pdf")
