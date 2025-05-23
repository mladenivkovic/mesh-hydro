/* Piecewise linear advection scheme */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/

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

  cell* c;  /* this cell */
  cell* uw; /* the cell upwind of the flux at the interface of *c */
  cell* dw; /* the cell downwind of the flux at the interface of *c */

#if NDIM == 1

  for (int i = BC; i < pars.nx + BC; i++) {
    c = &grid[i];
    if (c->prim.u[0] > 0) { /* we do upwind differencing */
      dw = c;
      uw = &(grid[i - 1]);
    } else {
      dw = &(grid[i + 1]);
      uw = c;
    }
    solver_compute_cell_pair_flux(c, uw, dw, dt, /*dimension=*/0);
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
        solver_compute_cell_pair_flux(c, uw, dw, dt, dimension);
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
        solver_compute_cell_pair_flux(c, uw, dw, dt, dimension);
      }
    }
  }
#endif /* ndim */
}

void solver_compute_cell_pair_flux(
  cell* c, cell* uw, cell* dw, const float* dt, int dim
) {
  /* --------------------------------------------------------------------
   * Compute the net flux for a given cell w.r.t. a specific cell pair
   * c:   pointer to cell to work with
   * uw:  the cell of the pair which is upwind
   * dw:  the cell of the cell pair which is downwind
   *
   * dim: integer along which dimension to advect. 0: x. 1: y.
   * -------------------------------------------------------------------- */

  pstate sl, sr; /* slopes left and right of the cell */

  gas_init_pstate(&sl);
  gas_init_pstate(&sr);

  limiter_get_advection_slope_left(c, &sl, dim);
  limiter_get_advection_slope_right(c, &sr, dim);

  /* assign the slopes to the upwind/downwind cells */
  pstate su; /* slope of upwind cell */
  pstate sd; /* slope of downwind cell */
  if (uw->id < dw->id) {
    su = sl;
    sd = sr;
  } else {
    su = sr;
    sd = sl;
  }

  /* Now compute the net fluxes */

  float dsu, dsd; /* dx - |v|*dt for upwind and downwind flux */

  if (uw->prim.u[dim] >= 0) {
    dsu = pars.dx - uw->prim.u[dim] * (*dt);
  } else {
    dsu = -pars.dx - uw->prim.u[dim] * (*dt);
  }
  if (dw->prim.u[dim] >= 0) {
    dsd = pars.dx - dw->prim.u[dim] * (*dt);
  } else {
    dsd = -pars.dx - dw->prim.u[dim] * (*dt);
  }

  c->pflux.rho = uw->prim.u[dim] * (uw->prim.rho + 0.5 * su.rho * dsu)
                 - dw->prim.u[dim] * (dw->prim.rho + 0.5 * sd.rho * dsd);
#ifndef ADVECTION_KEEP_VELOCITY_CONSTANT
  c->pflux.u[0] = uw->prim.u[dim] * (uw->prim.u[0] + 0.5 * su.u[0] * dsu)
                  - dw->prim.u[dim] * (dw->prim.u[0] + 0.5 * sd.u[0] * dsd);
  c->pflux.u[1] = uw->prim.u[dim] * (uw->prim.u[1] + 0.5 * su.u[1] * dsu)
                  - dw->prim.u[dim] * (dw->prim.u[1] + 0.5 * sd.u[1] * dsd);
#endif
  c->pflux.p = uw->prim.u[dim] * (uw->prim.p + 0.5 * su.p * dsu)
               - dw->prim.u[dim] * (dw->prim.p + 0.5 * sd.p * dsd);
}
