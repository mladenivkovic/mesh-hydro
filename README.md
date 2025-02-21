Mesh-based Hydro Solvers
===============================


Description
-----------------------

A program to learn about and play with mesh/cell/grid based hydrodynamics and
advection solvers. Currently, only finite volume methods are implemented. You
can pick between:

- Advection solvers:
    - piecewise constant advection
    - piecewise linear advection
    - weighted average flux (WAF) advection
- Hydro solvers:
    - Godunov (upwind)
    - Weighted Average Flux (WAF)
    - MUSCL-Hancock
- Riemann solvers:
    - Exact	(exact iterative solver)
    - HLLC	(Harten, Lax, van Leer approximate solver with contact wave)
    - TRRS	(Two Rarefaction approximate Riemann Solver)
    - TSRS	(Two Shock approximate Riemann Solver)
- Limiters:
    - no limiter
    - minmod
    - superbee
    - monotonized centered difference (MC)
    - van Leer
- Number of dimensions:
    - 1 or 2
- Source terms:
    - constant cartesian
    - constant radial w.r.t. box center
- Integrators (for source terms):
    - Runge Kutta 2
    - Runge Kutta 4

at compile time by setting the corresponding values in the Makefile.

To see what the code is able to do, you can for example run the test script in
`./program/test` and have a look at the resulting
`./program/test/test_results.pdf` (You will need to install the accompanying
`mesh_hydro_utils` python package for this to work. See the instructions below.)
Another option is to just run `./program/bin/run.sh` to run a Kelvin-Helmholtz
instability simulation.

You can read up how the code does that with the tex and pdf files in `tex/`.

The purpose of this project is to learn about, play with, and eventually teach
the basics of (finite volume) fluid dynamics. With that in mind, the project is
deliberately written to be easy to read, and well documented (check out
`/tex/equations/equations_and_implementation_details.pdf` !) So things are kept
simple and non-optimized.





Contents
----------------------------

- `./IC`: A collection of default initial condition files.
- `./program`: contains the actual software. The source code is in
  `./program/src`, the Makefile in `./program/bin`
- `./scripts` A collection of scripts; both for `python` and `bash`. Mainly
  intended to run a set of simulations in a coordinated fashion, produce
  comparisons etc. For plotting scripts, have a look at the scripts in the
  `mesh_hydro_utils` subdirectory/git submodule.
- `./python_module`: A git submodule containing the
  [mesh_hydro_utils](https://github.com/mladenivkovic/mesh_hydro_utils) python
  module. It contains convenience functions to generate initial conditions, plot
  ICs and outputs, and a Riemann solver. Note that you need to install it first
  for it to work. Instructions are given
  [below](#Getting-And-Installing-The-Python-Module).
- `./tex`: TeX documents containing every equations that is used, some
  implementation details, etc., some results and conclusions.



Installation
----------------------

### Requirements

- A good old C compiler. Code is written in C11 standard. I tested it with `gcc
  14.1.0` and `oneapi 2024`.
- GNU make to compile without much hassle. The Makefile is set up for GCC, but
  options for OneAPI is available.
- `python 3` with `numpy` and `matplotlib` for plotting.
- LaTeX to create the TeX files. I hardcoded the `pdflatex` command in the
  scripts. It doesn't require any fancy LaTeX packages.


### Getting The Code

You can get the code from the github repository:

```
$ git clone https://github.com/mladenivkovic/mesh-hydro.git
```

or

```
$ git clone git@github.com:mladenivkovic/mesh-hydro.git
```


### Getting And Installing The Python Module

The entire python module is stored within this repository as a git submodule of
its [own repository](https://github.com/mladenivkovic/mesh_hydro_utils).

Once you've cloned the mesh-hydro repository, you'll also need to tell git to
grab the submodules using

```
$ git submodule init
$ git submodule update
```

When completed successfully, the directory `./python_module` should now contain
some files. We now need to install this python module.

The easiest way is to navigate into the directory and install it locally using
e.g. `pip`:

```
$ cd python_module
$ pip install -e .
```

Alternatively (***albeit very discouraged***), you can add the directory
`./python_module/mesh_hydro_utils` to your `$PYTHONPATH`.







Usage
------------

```
cd program/bin/
make
./hydro paramfile ic-file
```

An example parameter file is given in `bin/example-paramfile.txt`.
Example IC files are in `IC/`


There is an option to just use the code as a Riemann solver.
To do that, use the `Makefile-Riemann` makefile in `/program/bin/`, i.e.

```
cd program/bin/
make -f Makefile-Riemann
./riemann paramfile ic-file
```

in the parameter file, you need to specify `nx` and `tmax` only, pretty much all
other parameters are ignored.








Things to keep in mind
----------------------------

- General:
    - You can only solve for one hyperbolic conservatio law at a time: You need to decide whether
      you're doing hydro or advection, and which solver you want to use.
    - All macros are in `defines.h`
    - Boxsize is always assumed to be 1 (in any dimension), and starts at zero. Can be changed in
      `defines.h` though.
    - There is only a uniform grid. If you want AMR or non-Cartesian geometries, do it yourself.
- Advection:
    - No advection solver uses a Riemann solver.
    - Piecewise constant advection can't apply limiters.
    - Reflective boundary conditions make no sense for any advection scheme. Transmissive ones
      aren't really useful either, so you probably should always go with periodic ones when dealing
      with advection.
- Hydrodynamics related:
    - This code is written assuming ideal gasses with adiabatic index gamma = 5/3. If you want
      something different, change it in `defines.h`
    - If you change the adiabatic index in the definitions, and want to overplot exact Riemann
      solutions using the python scripts, don't forget to change gamma in
      `py/module/hydro_riemann.py` as well!
    - WAF methods (advection and hydro) without the use of flux limiters introduces strong
      unphysical oscillations (as it is expected to). You probably won't be able to run the code
      over a long time interval. In fact, most vacuum examples can't even pass the tests.
    - The MUSCL-Hancock method doesn't work well with the TRRS solver. Like the WAF methods,
      oscillations occur if no limiter is used, and the code may crash and/or produce NANs.
    - The MUSCL-Hancock method doesn't handle vacuum very well with this interpretation. Best to
      avoid it.
    - There is no MC limiter implemented that works with the MUSCL-Hancock method. The code will
      throw an error if you try to run with it.
- Riemann related:
    - You can use the code as a Riemann solver only. To do that, use the `Makefile-Riemann` makefile
      in `/program/bin/`
    - ~~You can't use the HLLC solver as an individual solver to just solve a Riemann problem.~~ Now
      you can :)
    - The TRRS solver doesn't do well with the MUSCL method. Best to avoid it.
- Source terms:
    - Source terms have only been implemented for hydro applications. It should
      be straightforward to add them to advection though.
- Integrators:
    - Integrators are only employed if there are source terms to add to the
      Euler equations.








Parameter File
---------------------

The parameter file is written in plain text in the format
```
/* this is a comment and is ignored in the parameter file */
// this type of comment is ignored too

/* empty lines, like the one above and below this line, are ignored too */

<name> = <value>
```





Accepted names, their datatypes and default values are:


### Behaviour Options and Parameters

| name          |  default value    | type  | description                                                                   |
|---------------|-------------------|-------|-------------------------------------------------------------------------------|
| `verbose`     | = 0               | `int` | How talkative the code should be. 0 = quiet, 1 = talky, 2 = no secrets, 3 = debugging        |
|               |                   |       |                                                                               |
| `nstep_log`   | = 0               | `int` | Write log messages only ever `nstep_log` steps. If 0, will write every step.  |
|               |                   |       |                                                                               |



### Simulation Related Options and Parameters

| name              |  default value    | type  | description                                                                   |
|-------------------|-------------------|-------|-------------------------------------------------------------------------------|
| `nx`              | = 100             | `int` | Number of cells to use if you're running with a two-state type IC file. Otherwise, it needs to be specified in the initial conditions.  If you're not using a two-state IC, the value will be overwritten by the value given in the IC file.  |
|                   |                   |       |                                                                               |
| `ccfl`            | = 0.9             |`float`| courant factor; `dt = ccfl * dx / vmax`                                       |
|                   |                   |       |                                                                               |
| `nsteps`          | = 1               | `int` | Up to how many steps to do. If = 0, run until `t >= tmax`                     |
|                   |                   |       |                                                                               |
| `tmax`            | = 0               |`float`|  Up to which time to simulate. If `nsteps` is given, will stop running if `nsteps` steps are reached before `tmax` is.     |
|                   |                   |       |                                                                               |
| `force_dt`        | = 0               |`float`| force a time step size. If a smaller time step is required, the sim will stop.|
|                   |                   |       |                                                                               |
| `boundary`        | = 0               | `int` | Boundary conditions  0: periodic. 1: reflective. 2: transmissive. This sets the boundary conditions for all walls. |
|                   |                   |       |                                                                               |



### Output related Options and Parameters


| name          |  default value    | type    | description                                                                   |
|---------------|-------------------|---------|-------------------------------------------------------------------------------|
| `foutput`     | = 0               | `int`   | Frequency of writing outputs in number of steps. If = 0, will only write initial and final steps.  |
|               |                   |         |                                                                               |
| `dt_out`      | = 0               |`float`  | Frequency of writing outputs in time intervals. Code will always write initial and final steps as well.  |
|               |                   |         |                                                                               |
| `toutfile`    | None              |`string` | File name containing desired times (in code units) of output. Syntax of the file: One float per line with increasing value.  |
|               |                   |         |                                                                               |
| `basename`    | None              |`string` | Basename for outputs.  If not given, a basename will be generated based on compilation parameters and IC filename.       |
|               |                   |         |                                                                               |



### Source Term related Options and Parameters

The source term related options will only take effect if the code has been compiled to add source terms.


| name              |  default value    | type   | description                                                                   |
|-------------------|-------------------|--------|-------------------------------------------------------------------------------|
| `src_const_acc_x` | = 0               | `float`| constant acceleration in x direction for constant source terms                |
|                   |                   |        |                                                                               |
| `src_const_acc_y` | = 0               | `float`| constant acceleration in y direction for constant source terms                |
|                   |                   |        |                                                                               |
| `src_const_acc_r` | = 0               | `float`| constant acceleration in radial direction for radial source terms             |
|                   |                   |        |                                                                               |


Initial Conditions
------------------------------------

- The program reads two types of IC files.
- In any case, they're expected to be formatted text.
- In both IC file types, lines starting with `//` or `/*` will be recognized as
  comments and skipped. Empty lines are skipped as well.
- Some example python scripts that generate initial conditions are given in
  `./python_module/scripts/IC`. Note that for the directory `./python_module/`
  to contain any files, you first need to initialize the submodule. See the
  instructions [above](#Getting-And-Installing-The-Python-Module).




### Two-state ICs

You can use a Riemann-problem two-state initial condition file as follows:

```
filetype = two-state
rho_L   = <float>
u_L     = <float>
p_L     = <float>
rho_R   = <float>
u_R     = <float>
p_R     = <float>
```

The line
```
filetype = two-state
```

**must** be the first non-comment non-empty line. The order of the other
parameters is arbitrary, but they must be named `rho_L`, `u_L`, `p_L`, `rho_R`,
`u_R`, `p_R`.

The discontinuity between the changes will be in the middle along the x axis.
The coordinates will be printed to screen.

If the code is supposed to run in 2D, then the split will be along the x axis as
well, and just copied along the y axis. Fluid velocity in y direction will be
set to zero, `u_L` and `u_R` will be set as `u_x`.



### Arbitrary ICs

You can provide an individual value for density, velocity, and pressure for each
cell. The IC file format is:

The lines
```
filetype = arbitrary
nx = <int>
ndim = <int>
```

**must** be the first non-comment non-empty lines, in that order.

**FOR 1D:**

```
filetype = arbitrary
nx = <integer, number of cells in any dimension>
ndim = 1
<density in cell 0> <velocity in cell 0> <pressure in cell 0>
<density in cell 1> <velocity in cell 1> <pressure in cell 1>
                            .
                            .
                            .
<density in cell nx-1> <velocity cell nx-1> <pressure in cell nx-1>
```

`cell 0` is the leftmost cell. All values for density, velocity, and pressure
must be floats. You can put comments and empy lines wherever you feel like it.






**FOR 2D:**

```
filetype = arbitrary
nx = <integer, number of cells in any dimension>
ndim = 2
<density in cell (0, 0)> <x velocity in cell (0, 0)> <y velocity in cell (0, 0)> <pressure in cell (0, 0)>
<density in cell (1, 0)> <x velocity in cell (1, 0)> <y velocity in cell (1, 0)> <pressure in cell (1, 0)>
                                     .
                                     .
                                     .
<density in cell (nx-1, 0)> <x velocity cell (nx-1, 0)> <y velocity in cell (nx-1, 0)> <pressure in cell (nx-1, 0)>
<density in cell (0, 1)> <x velocity in cell (0, 1)> <y velocity in cell (0, 1)> <pressure in cell (0, 1)>
<density in cell (1, 1)> <x velocity in cell (1, 1)> <y velocity in cell (1, 1)> <pressure in cell (1, 1)>
                                     .
                                     .
                                     .
<density in cell (nx-1, 1)> <x velocity cell (nx-1, 1)> <y velocity in cell (nx-1, 1)> <pressure in cell (nx-1, nx-1)>
                                     .
                                     .
                                     .
<density in cell (0, nx-1)> <x velocity in cell (0, nx-1)> <y velocity in cell (0, nx-1)> <pressure in cell (0, nx-1)>
<density in cell (1, nx-1)> <x velocity in cell (1, nx-1)> <y velocity in cell (1, nx-1)> <pressure in cell (1, nx-1)>
                                     .
                                     .
                                     .
<density in cell (nx-1, nx-1)> <x velocity cell (nx-1, nx-1)> <y velocity in cell (nx-1, nx-1)> <pressure in cell (nx-1, nx-1)>
```

`cell (0, 0)` is the lower left corner of the box. First index is x direction, second is y. All values for density, velocity, and pressure must be floats.
You can put comments and empy lines wherever you feel like it.










Output
---------------------


If no `basename` is given in the parameter file, the output file name will be
generated as follows:

`<ICfile-without-suffix>-<SOLVER>-<RIEMANN-SOLVER>-<LIMITER>-<NDIM>D-<snapshot nr>.out`

e.g.

`run-ADVECTION-NO_LIMITER-2D-0001.out`

The output files are written in formatted text, and their content should be
self-explainatory:




**FOR 1D**:

```
# ndim =  1
# nx =    <number of cells used>
# t =     <current time, float>
# nsteps =  <current step of the simulation>
#            x          rho            u            p
< x value of cell 0 > <density of cell 0> <velocity of cell 0> < pressure of cell 0>
                            .
                            .
                            .
                            .
< x value of cell nx-1 > <density of cell nx-1> <velocity of cell nx-1> < pressure of cell nx-1>
```





**FOR 2D:**

```
# ndim =  2
# nx =    <number of cells used>
# t =     <current time, float>
# nsteps =  <current step of the simulation>
#            x            y          rho          u_x          u_y            p
<x value of cell (0, 0)> <y value of cell (0, 0)> <density in cell (0, 0)> <x velocity in cell (0, 0)> <y velocity in cell (0, 0)> <pressure in cell (0, 0)>
<x value of cell (1, 0)> <y value of cell (1, 0)> <density in cell (1, 0)> <x velocity in cell (1, 0)> <y velocity in cell (1, 0)> <pressure in cell (1, 0)>
                                                 .
                                                 .
                                                 .
<x value of cell (nx-1, 0)> <y value of cell (nx-1, 0)> <density in cell (nx-1, 0)> <x velocity cell (nx-1, 0)> <y velocity in cell (nx-1, 0)> <pressure in cell (nx-1, 0)>
<x value of cell (0, 1)> <y value of cell (0, 1)> <density in cell (0, 1)> <x velocity in cell (0, 1)> <y velocity in cell (0, 1)> <pressure in cell (0, 1)>
<x value of cell (1, 1)> <y value of cell (1, 1)> <density in cell (1, 1)> <x velocity in cell (1, 1)> <y velocity in cell (1, 1)> <pressure in cell (1, 1)>
                                                 .
                                                 .
                                                 .
<x value of cell (nx-1, 1)> <y value of cell (nx-1, 1)> <density in cell (nx-1, 1)> <x velocity cell (nx-1, 1)> <y velocity in cell (nx-1, 1)> <pressure in cell (nx-1, nx-1)>
                                                 .
                                                 .
                                                 .
<x value of cell (0, nx-1)> <y value of cell (0, nx-1)> <density in cell (0, nx-1)> <x velocity in cell (0, nx-1)> <y velocity in cell (0, nx-1)> <pressure in cell (0, nx-1)>
<x value of cell (1, nx-1)> <y value of cell (1, nx-1)> <density in cell (1, nx-1)> <x velocity in cell (1, nx-1)> <y velocity in cell (1, nx-1)> <pressure in cell (1, nx-1)>
                                                 .
                                                 .
                                                 .
<x value of cell (nx-1, nx-1)> <y value of cell (nx-1, nx-1)> <density in cell (nx-1, nx-1)> <x velocity cell (nx-1, nx-1)> <y velocity in cell (nx-1, nx-1)> <pressure in cell (nx-1, nx-1)>
```







Visualisation
-----------------------------------

Some basic scripts to visualize ICs and outputs are given in the
`./python_module/scripts/plotting` directory. See the `README.md` in the
`./python_module/scripts` directory for more details. Note that for the
directory `./python_module/` to contain any files, you first need to initialize
the submodule. See the instructions [above](#Getting-And-Installing-The-Python-Module).




Tinkering with the Code
----------------------------------

I tried keeping the code as modular as possible, so adding/removing stuff should
work fine. If you want to do that, here's a list of functions that could be
useful:

### in `utils.h`:

- `log_extra()`, `debugmessage()`, `throw_error()`: Functions that will print to
  stdout depending on the level of verbosity you set in the paramfile.
- `printbool()`: print boolean as true or false to stdout.




### in `cell.h`:

- `cell_print_grid()`: prints out chosen grid quantity for the entire grid to
  stdout. Works independently of dimension of the code, so you can always call
  it.
- `cell_print_grid_part()`: prints out chosen grid quantity for chosen cell
  intervals, i.e. give it xmin and xmax etc.






Contributing
-----------------------------------

Yes please!
