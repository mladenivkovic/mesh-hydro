/* Piecewise constant advection scheme */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/

#include <math.h>
#include <stdio.h>

#include "cell.h"
#include "defines.h"
#include "io.h"
#include "params.h"
#include "solver.h"
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
  solver_get_advection_dt(dt);
  /* check this here in case you need to limit time step for output */
  *write_output = io_is_output_step(*t, dt, step);

#if NDIM == 1

  solver_compute_fluxes(/*dimension =*/0);
  solver_advance_step_advection(dt);

#elif NDIM == 2

  int dimension = step % 2; /* gives 0 or 1, switching each step */
  solver_compute_fluxes(dimension);
  solver_advance_step_advection(dt);

  dimension = (dimension + 1) % 2; /* 1 -> 0 or 0 -> 1 */
  solver_init_step();
  solver_compute_fluxes(dimension);
  solver_advance_step_advection(dt);

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
  cell_set_boundary();
}

void solver_compute_fluxes(int dimension) {
  /* ------------------------------------------------------
   * Computes the actual fluxes between cells
   * Here we compute F^{n+1/2}_{i-1/2} - F^{n+1/2}_{i+1/2}
   * and store it in cell.pflux
   * int dimension: 0 for x, 1 for y. Only used when strang
   * splitting is employed.
   * ------------------------------------------------------ */

  debugmessage("Called solver_compute_fluxes; dimension = %d", dimension);

  cell *c;  /* this cell */
  cell *uw; /* the cell upwind of the flux at the interface of *c */
  cell *dw; /* the cell downwind of the flux at the interface of *c */

#if NDIM == 1

  for (int i = BC; i < pars.nx + BC; i++) {
    c = &(grid[i]);
    if (c->prim.u[0] > 0) { /* we do upwind differencing */
      dw = c;
      uw = &(grid[i - 1]);
    } else {
      dw = &(grid[i + 1]);
      uw = c;
    }
    solver_compute_cell_pair_flux(c, uw, dw, /*dimension=*/0);
  }

#elif NDIM == 2

  if (dimension == 0) {
    for (int i = BC; i < pars.nx + BC; i++) {
      for (int j = BC; j < pars.nx + BC; j++) {
        c = &(grid[i][j]);
        if (c->prim.u[0] > 0) {
          dw = c;
          uw = &(grid[i - 1][j]);
        } else {
          dw = &(grid[i + 1][j]);
          uw = c;
        }
        solver_compute_cell_pair_flux(c, uw, dw, dimension);
      }
    }
  } else if (dimension == 1) {
    for (int i = BC; i < pars.nx + BC; i++) {
      for (int j = BC; j < pars.nx + BC; j++) {
        c = &(grid[i][j]);
        if (c->prim.u[1] > 0) {
          dw = c;
          uw = &(grid[i][j - 1]);
        } else {
          dw = &(grid[i][j + 1]);
          uw = c;
        }
        solver_compute_cell_pair_flux(c, uw, dw, dimension);
      }
    }
  }
#endif /* ndim */
}

void solver_compute_cell_pair_flux(cell *c, cell *uw, cell *dw, int dim) {
  /* --------------------------------------------------------------------
   * Compute the net flux for a given cell w.r.t. a specific cell pair
   * c:   pointer to cell to work with
   * uw:  the cell of the pair which is upwind
   * dw:  the cell of the cell pair which is downwind
   *
   * dim: integer along which dimension to advect. 0: x. 1: y.
   * -------------------------------------------------------------------- */

  c->pflux.rho =
      uw->prim.rho * uw->prim.u[dim] - dw->prim.rho * dw->prim.u[dim];
#ifndef ADVECTION_KEEP_VELOCITY_CONSTANT
  c->pflux.u[0] =
      uw->prim.u[0] * uw->prim.u[dim] - dw->prim.u[0] * dw->prim.u[dim];
  c->pflux.u[1] =
      uw->prim.u[1] * uw->prim.u[dim] - dw->prim.u[1] * dw->prim.u[dim];
#endif
  c->pflux.p = uw->prim.p * uw->prim.u[dim] - dw->prim.p * dw->prim.u[dim];
}
