/* First order Godunov scheme */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/

#include <math.h>
#include <stdio.h>

#include "cell.h"
#include "defines.h"
#include "io.h"
#include "limiter.h"
#include "params.h"
#include "riemann.h"
#include "solver.h"
#include "sources.h"
#include "utils.h"

#if NDIM == 1
extern cell *grid;
#elif NDIM == 2
extern cell **grid;
#endif

extern params pars;

void solver_step(float *t, float *dt, int step, int *write_output) {
  /* -------------------------------------------------------
   * Main routine for the actual hydro step
   * ------------------------------------------------------- */

  solver_init_step();
  solver_get_hydro_dt(dt, step);
  /* check this here in case you need to limit time step for output */
  *write_output = io_is_output_step(*t, dt, step);

#if NDIM == 1

  solver_compute_fluxes(dt, /*dimension =*/0);
  solver_advance_step_hydro(dt, /*dimension=*/0);

#elif NDIM == 2

  int dimension = step % 2; /* gives 0 or 1, switching each step */
  solver_compute_fluxes(dt, dimension);
  solver_advance_step_hydro(dt, dimension);

  dimension = (dimension + 1) % 2; /* 1 -> 0 or 0 -> 1 */
  solver_init_step();
  solver_compute_fluxes(dt, dimension);
  solver_advance_step_hydro(dt, dimension);
  /* cell_reset_fluxes(); */ /* will be done in solver_init_step */

#endif /* ndim = 2 */

  /* after the step is finished, compute the source terms if there are
   * any present. This gives a first order accurate scheme */
#ifdef WITH_SOURCES
  sources_update_state(*dt);
#endif
}

void solver_init_step() {
  /* ---------------------------------------------
   * Do everything that needs to be done before
   * we can compute the fluxes, the timestep, and
   * finally advance the simulation
   * --------------------------------------------- */

  debugmessage("Called solver_init_step");
  cell_reset_fluxes();
  cell_get_pstates_from_cstates();
  cell_set_boundary();
}

void solver_compute_fluxes(float *dt, int dimension) {
  /* ------------------------------------------------------------
   * Compute the flux F_{i+1/2} (or G_{i+1/2} if dimension == 1)
   * and store it in cell i.
   * BUT: Start with the first boundary cell, such that even the
   * first real cell can access F_{i-1/2} by accessing the
   * neighbour at i-1
   * ----------------------------------------------------------- */

  debugmessage("Called solver_compute_fluxes; dimension = %d", dimension);

  cell *left;  /* this cell */
  cell *right; /* the right neighbour */

#if NDIM == 1

  for (int i = BC - 1; i < pars.nx + BC; i++) {
    left = &grid[i];
    right = &grid[i + 1];
    solver_compute_cell_pair_flux(left, right, dt, /*dimension=*/0);
  }

#elif NDIM == 2

  if (dimension == 0) {
    for (int i = BC - 1; i < pars.nx + BC; i++) {
      for (int j = BC - 1; j < pars.nx + BC; j++) {
        left = &(grid[i][j]);
        right = &(grid[i + 1][j]);
        /* debugmessage("Calling solver_compute_cell_pair_flux for cell %d %d",
         * i, j); */
        solver_compute_cell_pair_flux(left, right, dt, dimension);
      }
    }
  } else if (dimension == 1) {
    for (int i = BC - 1; i < pars.nx + BC; i++) {
      for (int j = BC - 1; j < pars.nx + BC; j++) {
        left = &(grid[i][j]);
        right = &(grid[i][j + 1]);
        /* debugmessage("Calling solver_compute_cell_pair_flux for cell %d %d",
         * i, j); */
        solver_compute_cell_pair_flux(left, right, dt, dimension);
      }
    }
  }

#endif /* ndim */
}

void solver_compute_cell_pair_flux(cell *left, cell *right, float *dt,
                                   int dim) {
  /* --------------------------------------------------------------------
   * Compute the net flux for a given cell w.r.t. a specific cell pair
   * left:  pointer to cell which stores the left state
   * right: pointer to cell which stores the right state
   * dt:    current time step
   * dim:   integer along which dimension to advect. 0: x. 1: y.
   *
   * Here, we just solve the Riemann problem between the left and right
   * primitive states, and sample the solution at x = x/t = 0, because
   * that is where we set the initial boundary in the local coordinate
   * system between the left and right cell.
   * -------------------------------------------------------------------- */

#if RIEMANN == HLLC
  /* the HLLC solver gives us the flux directly. */
  riemann_solve_hllc(&left->prim, &right->prim, &left->cflux, /*xovert=*/0.0,
                     dim);
#else
  pstate solution;
  gas_init_pstate(&solution);

  /* solve the Riemann problem */
  riemann_solve(&left->prim, &right->prim, &solution, /*xovert=*/0.0, dim);

  /* from the primitive states, compute and store F_{i+1/2} */
  gas_get_cflux_from_pstate(&solution, &left->cflux, dim);
#endif
}
