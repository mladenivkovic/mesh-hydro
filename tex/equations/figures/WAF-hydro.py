#!/usr/bin/env python3

#-----------------------------------------------------------------
# Create a plot for the hydro WAF method
#-----------------------------------------------------------------


import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.patches as patches

params = {
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'font.size': 14,
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

xmin = -10
xmax = -xmin
dx = (xmax - xmin) * 0.01
ymin = 0
ymax = 10

fig = plt.figure(figsize=(6, 3))
ax = fig.add_subplot(111)
ax.axis("off")

ax.plot([xmin-dx, xmax+dx], [0, 0], 'k', zorder=10) # x axis
ax.plot([0, 0], [ymin, ymax], 'k', zorder=10)       # t axis

ax.plot([xmin-dx, xmax+dx], [4.8, 4.8], 'k', ls=':', lw=1.) # dt/2 line
ax.plot([xmin-dx, xmax+dx], [9.6, 9.6], 'k', ls=':', lw=1.) # dt line
ax.plot([xmax/2, xmax/2], [ymin, ymax+dx], 'k', ls=':', lw=1.) # delta x line
ax.plot([-xmax/2, -xmax/2], [ymin, ymax+dx], 'k', ls=':', lw=1.) # - delta x line




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


# characteristics
ax.plot([0, x_for_ymax(-1.5)], [0, ymax], 'k', lw=1)
ax.plot([0, x_for_ymax(5)], [0, ymax], 'k', lw=1)
ax.plot([0, x_for_ymax(1.5)], [0, ymax], 'k', lw=1)


# point A0
ax.scatter([xmin/2], [4.8], c='k', s=20, zorder=1e6)
plt.text(xmin/2-0.2, 4.8-0.2, r"$A_0$",
    verticalalignment="top",
    horizontalalignment='right',)

# point A1
ax.scatter([x_for_y(4.8, -1.5)], [4.8], c='k', s=20, zorder=1e6)
plt.text(x_for_y(4.8, -1.5)-0.2, 4.8-0.2, "$A_1$",
    color='k',
    verticalalignment="top",
    horizontalalignment='right',)

# point A2
ax.scatter([x_for_y(4.8, 5)], [4.8], c='k', s=20, zorder=1e6)
plt.text(x_for_y(4.8, 5)+0.2, 4.8-0.2, "$A_2$",
    color='k',
    verticalalignment="top",
    horizontalalignment='left',)

# point A3
ax.scatter([x_for_y(4.8, 1.5)], [4.8], c='k', s=20, zorder=1e6)
plt.text(x_for_y(4.8, 1.5)+0.2, 4.8-0.2, "$A_3$",
    color='k',
    verticalalignment="top",
    horizontalalignment='left',)

# point A4
ax.scatter([xmax/2], [4.8], c='k', s=20, zorder=1e6)
plt.text(xmax/2+0.2, 4.8-0.2, "$A_4$",
    verticalalignment="top",
    horizontalalignment='left',)

# U_L
plt.text(-7, 7, "$\mathbf{U}_L$",
    color='k',
    verticalalignment="top",
    horizontalalignment='left',)

# Ustar_L
plt.text(-3, 7, "$\mathbf{U}_L^*$",
    color='k',
    verticalalignment="top",
    horizontalalignment='left',)

# Ustar_R
plt.text(2, 7, "$\mathbf{U}_R^*$",
    color='k',
    verticalalignment="top",
    horizontalalignment='left',)

# U_R
plt.text(7, 7, "$\mathbf{U}_R$",
    color='k',
    verticalalignment="top",
    horizontalalignment='left',)



# S1
plt.text(-7, 11, "$S_1$",
    color='k',
    verticalalignment="top",
    horizontalalignment='left',)

# S2
plt.text(2, 11, "$S_2$",
    color='k',
    verticalalignment="top",
    horizontalalignment='left',)

# S3
plt.text(7, 11, "$S_3$",
    color='k',
    verticalalignment="top",
    horizontalalignment='left',)

# origin
plt.text(0, -1, r"$x_{i+1/2}$",
    verticalalignment="center",
    horizontalalignment='left',)
# x label
plt.text(10.5, 0.0, r"$x$",
    verticalalignment="center",
    horizontalalignment='left',)
# t label
plt.text(0, 11, r"t",
    verticalalignment="center",
    horizontalalignment='center',)
# Dt/2
plt.text(10.5, 4.8, r"$\Delta t$/2",
    verticalalignment="center",
    horizontalalignment='left',)
# Dt
plt.text(10.5, 9.6, r"$\Delta t$",
    verticalalignment="center",
    horizontalalignment='left',)
# -Dx/2
plt.text(xmin/2, -1, r"$-\Delta x$/2",
    verticalalignment="center",
    horizontalalignment='center',)
# Dx/2
plt.text(xmax/2, -1, r"$\Delta x$/2",
    verticalalignment="center",
    horizontalalignment='center',)



# U_i coloring
ax.add_patch(patches.Rectangle((xmin, ymin), width=(xmax-xmin)/2, height=(ymax-ymin), fill=True, facecolor='lightblue'))
# U_i+1 coloring
ax.add_patch(patches.Rectangle((0, ymin), width=(xmax-xmin)/2, height=(ymax-ymin), fill=True, facecolor='thistle'))


#  plt.show()
plt.savefig("WAF-hydro.pdf")
