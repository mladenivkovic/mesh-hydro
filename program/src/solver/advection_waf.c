/* Weighted Average Flux advection scheme

 * Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/

#include <math.h>

#include "cell.h"
#include "defines.h"
#include "io.h"
#include "limiter.h"
#include "params.h"
#include "solver.h"
#include "utils.h"

extern params pars;

void solver_step(const float* t, float* dt, int step, int* write_output) {
  /* -------------------------------------------------------
   * Main routine for the actual hydro step
   * ------------------------------------------------------- */

  solver_init_step();
  solver_get_advection_dt(dt);
  /* check this here in case you need to limit time step for output */
  *write_output = io_is_output_step(*t, dt, step);

#if NDIM == 1

  solver_compute_fluxes(dt, /*dimension =*/0);
  solver_advance_step_advection(dt);

#elif NDIM == 2

  int dimension = step % 2; /* gives 0 or 1, switching each step */
  solver_compute_fluxes(dt, dimension);
  solver_advance_step_advection(dt);

  dimension = (dimension + 1) % 2; /* 1 -> 0 or 0 -> 1 */
  solver_init_step();
  solver_compute_fluxes(dt, dimension);
  solver_advance_step_advection(dt);

#endif
}

void solver_init_step(void) {
  /* ---------------------------------------------
   * Do everything that needs to be done before
   * we can compute the fluxes, the timestep, and
   * finally advance the simulation
   * --------------------------------------------- */

  debugmessage("Called solver_init_step");
  cell_reset_fluxes();
  cell_set_boundary();
}

void solver_compute_fluxes(float* dt, int dimension) {
  /* ------------------------------------------------------
   * Computes the actual *net* fluxes between cells
   * Here we compute F^{n+1/2}_{i-1/2} - F^{n+1/2}_{i+1/2}
   * and store it in cell.pflux
   * dimension: 0 for x, 1 for y. Only used with Strang
   * splitting.
   * ------------------------------------------------------ */

  debugmessage("Called solver_compute_fluxes; dimension = %d", dimension);

  cell* c; /* this cell */
  cell* n; /* the neighbour we're dealing with, i.e. i+1 or j+1 */

#if NDIM == 1

  /* computes F_{i+1/2}, so we need to start at index BC-1:
   * F_{i+1/2} = F_{(i+1) - 1/2}*/
  for (int i = BC - 1; i < pars.nx + BC; i++) {
    c = &grid[i];
    n = &grid[i + 1];
    solver_compute_cell_pair_flux(c, n, dt, /*dim=*/0);
  }

#elif NDIM == 2

  if (dimension == 0) {
    /* computes F_{i+1/2}, so we need to start at index BC-1:
     * F_{i+1/2} = F_{(i+1) - 1/2}*/
    for (int i = BC - 1; i < pars.nx + BC; i++) {
      for (int j = BC - 1; j < pars.nx + BC; j++) {
        c = &(grid[i][j]);
        n = &(grid[i + 1][j]);
        solver_compute_cell_pair_flux(c, n, dt, dimension);
      }
    }
  } else if (dimension == 1) {
    /* computes F_{i+1/2}, so we need to start at index BC-1:
     * F_{i+1/2} = F_{(i+1) - 1/2}*/
    for (int i = BC - 1; i < pars.nx + BC; i++) {
      for (int j = BC - 1; j < pars.nx + BC; j++) {
        c = &(grid[i][j]);
        n = &(grid[i][j + 1]);
        solver_compute_cell_pair_flux(c, n, dt, dimension);
      }
    }
  }
#endif /* ndim */
}

void solver_compute_cell_pair_flux(cell* c, cell* n, const float* dt, int dim) {
  /* --------------------------------------------------------------------
   * Compute F_{i+1/2} where cell* c is cell with index i.
   * Then subtract the flux from the net flux of cell c, and add it to the
   * net flux of cell i+1.
   *
   * c:   pointer to cell to work with
   * n:   pointer to neighbour cell; i.e. cell with index i+1 or j+1
   *
   * dim: integer along which dimension to advect. 0: x. 1: y.
   * -------------------------------------------------------------------- */

  pstate phi;
  gas_init_pstate(&phi);
  limiter_get_phi(c, &phi, dim);

  float vel    = c->prim.u[dim];
  float abscfl = (*dt) / pars.dx * fabsf(vel);

  pstate psi;
  gas_init_pstate(&psi);
  psi.rho  = 1. - (1. - abscfl) * phi.rho;
  psi.u[0] = 1. - (1. - abscfl) * phi.u[0];
  psi.u[1] = 1. - (1. - abscfl) * phi.u[1];
  psi.p    = 1. - (1. - abscfl) * phi.p;

  float s = 1.; /* sign(velocity) */
  if (vel <= 0) { s = -1.; }

  float flux = 0.;

  flux = 0.5 * (1. + s * psi.rho) * vel * c->prim.rho
         + 0.5 * (1. - s * psi.rho) * vel * n->prim.rho;
  c->pflux.rho -= flux;
  n->pflux.rho += flux;
#ifndef ADVECTION_KEEP_VELOCITY_CONSTANT
  flux = 0.5 * (1. + s * psi.u[0]) * vel * c->prim.u[0]
         + 0.5 * (1. - s * psi.u[0]) * vel * n->prim.u[0];
  c->pflux.u[0] -= flux;
  n->pflux.u[0] += flux;
  flux = 0.5 * (1. + s * psi.u[1]) * vel * c->prim.u[1]
         + 0.5 * (1. - s * psi.u[1]) * vel * n->prim.u[1];
  c->pflux.u[1] -= flux;
  n->pflux.u[1] += flux;
#endif
  flux = 0.5 * (1. + s * psi.p) * vel * c->prim.p
         + 0.5 * (1. - s * psi.p) * vel * n->prim.p;
  c->pflux.p -= flux;
  n->pflux.p += flux;
}
