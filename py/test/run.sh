#!/bin/bash


# Test your scripts stored in ../plotting on ICs and outputs

# ../plotting/plot_all_density_individually.py advection-00*out
# ../plotting/plot_all_density_individually.py advection-2D-00*out
# ../plotting/plot_all_density.py advection-00*out
# ../plotting/plot_all_results.py advection-00*out
# ../plotting/plot_all_results.py advection-2D-00*out
# ../plotting/plot_all_results_individually.py advection-00*out
# ../plotting/plot_all_results_individually.py advection-2D-00*out
# ../plotting/plot_density.py advection-0009.out
# ../plotting/plot_density.py advection-2D-0004.out
# ../plotting/plot_IC_3D.py ic-1D.dat
# ../plotting/plot_IC_3D.py ic-2D.dat
# ../plotting/plot_IC_3D.py ic-twostate.dat
# ../plotting/plot_IC.py ic-1D.dat
# ../plotting/plot_IC.py ic-2D.dat
# ../plotting/plot_IC.py ic-twostate.dat
# ../plotting/plot_result_3D.py advection-0010.out
# ../plotting/plot_result_3D.py advection-2D-0004.out
# ../plotting/plot_result.py advection-0010.out
# ../plotting/plot_result.py advection-2D-0004.out
# ../plotting/plot_riemann_result.py sod-shock-0001.out ic-twostate.dat
# ../plotting/plot_all_riemann_results.py sod-shock-000*.out ic-twostate.dat
../plotting/solve_riemann.py  ic-twostate.dat 0.25

# ../plotting/artsy_plot.py advection-2D-0004.out
# ../plotting/movie_density.py advection-2D-0004.out
# ../plotting/plot_all_2D_velnorm_individually.py advection-2D-0004.out
