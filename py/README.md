Python Scripts
===========================


Just add the `./module` directory to your `PYTHONPATH` first!


Plotting
------------------

In the `./plotting` directory. In the following, "output" refers to output files created by running the hydro code.

- `movie_density.py`: Plot the density of all given output files. You need to fix some parameters manually in the script so that the images are standardized to make a movie out of it.
- `plot_all_density_individually.py`: Plots all results given in individual images. Only plots the density.
- `plot_all_density.py`: Overplots all results in one plot. Only works for 1D. Only plots the density.
- `plot_all_results_individually.py`: Plots all results given in individual images.
- `plot_all_results.py`: Overplots all results in one plot. Only works for 1D.
- `plot_density.py`: Plots output, density only.
- `plot_IC_3D.py`: Plot some initial conditions file. Plot 2D IC in a 3D surface plot.
- `plot_IC.py`: Plot some initial conditions file. It figures out what type of IC it is and what dimension it is by itself.
- `plot_result_3D.py`: Plot output written by the hydro code. Plot 2D results in a 3D surface plot.
- `plot_result.py`: Plot output written by the hydro code. It'll figure out the dimension etc by itself.
- `plot_riemann_result.py`: Plot the output of a two-state/Riemann problem, and overplot analytical solution.
- `solve_riemann.py`: Solve a riemann problem from a given IC at a given time and plot it.




Generating IC
------------------

In the `./IC` directory. Run them to generate some easy to do profiles etc.
`output_to_IC.py` generates an IC file from a result file.



Modules
------------------

In the `./module` directory.

Not meant to be run, but contain functions that other scripts call.

- `deprecated.py`:        Stuff that shouldn't be used any more.
- `hydro_io.py`:          reading/writing routines.
- `hydro_plotting.py`:    Plotting routines.
- `hydro_riemann.py`:     An exact Riemann solver.
- `hydro_utils.py`:       misc small utilities.



Evaluations
------------------

Contains plotting, comparing, and evaluation scripts for different tasks.
E.g. compute and plot the convergence of a class of methods depending on grid spacing, dt, ...
