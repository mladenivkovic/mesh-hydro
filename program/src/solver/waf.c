/* Weighted Average Flux Hydro Scheme  */

/* Written by Mladen Ivkovic, APR 2020
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

extern params pars;

void solver_step(const float *t, float *dt, int step, int *write_output) {
  /* -------------------------------------------------------
   * Main routine for the actual hydro step
   * ------------------------------------------------------- */

  solver_init_step();
  solver_get_hydro_dt(dt, step);
  /* check this here in case you need to limit time step for output */
  *write_output = io_is_output_step(*t, dt, step);

#if defined WITH_SOURCES || NDIM > 1
  float dthalf = 0.5 * (*dt);
#endif

#ifdef WITH_SOURCES
  /* if we have source terms to add to the Euler equations, update them
   * for half the timestep now, and for half the time step after the
   * homogeneous part of the Euler equations has been solved. This results
   * in a second order accurate scheme to include source terms.*/
  debugmessage("Computing source term ODE for half timestep before hydro step");
  sources_update_state(dthalf);
#endif

#if NDIM == 1

  solver_compute_fluxes(dt, /*dimension =*/0);
  solver_advance_step_hydro(dt, /*dimension=*/0);

#elif NDIM == 2

  int dim;

  dim = step % 2;
  debugmessage("Advancing dim=%d for half timestep", dim);
  solver_compute_fluxes(&dthalf, dim);
  solver_advance_step_hydro(&dthalf, dim);

  dim = (step + 1) % 2;
  debugmessage("Advancing dim=%d for full timestep", dim);
  solver_init_step();
  solver_compute_fluxes(dt, dim);
  solver_advance_step_hydro(dt, dim);

  dim = step % 2;
  debugmessage("Advancing dim=%d for half timestep", dim);
  solver_init_step();
  solver_compute_fluxes(&dthalf, dim);
  solver_advance_step_hydro(&dthalf, dim);

#endif

#ifdef WITH_SOURCES
  /* Update the source terms for the second time over half the time step
   * for second order accuracy. */
  debugmessage("Computing source term ODE for half timestep after hydro step "
               "has finished");
  sources_update_state(dthalf);
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
    solver_prepare_flux_computation(left, right, /*dimension=*/0);
  }

  for (int i = BC - 1; i < pars.nx + BC; i++) {
    left = &grid[i];
    right = &grid[i + 1];
    solver_compute_cell_pair_flux(left, right, dt, /*dim=*/0);
  }

#elif NDIM == 2

  if (dimension == 0) {
    for (int i = BC - 1; i < pars.nx + BC; i++) {
      for (int j = BC - 1; j < pars.nx + BC; j++) {
        left = &(grid[i][j]);
        right = &(grid[i + 1][j]);
        solver_prepare_flux_computation(left, right, dimension);
      }
    }
    for (int i = BC - 1; i < pars.nx + BC; i++) {
      for (int j = BC - 1; j < pars.nx + BC; j++) {
        left = &(grid[i][j]);
        right = &(grid[i + 1][j]);
        solver_compute_cell_pair_flux(left, right, dt, dimension);
      }
    }
  } else if (dimension == 1) {
    for (int i = BC - 1; i < pars.nx + BC; i++) {
      for (int j = BC - 1; j < pars.nx + BC; j++) {
        left = &(grid[i][j]);
        right = &(grid[i][j + 1]);
        solver_prepare_flux_computation(left, right, dimension);
      }
    }
    for (int i = BC - 1; i < pars.nx + BC; i++) {
      for (int j = BC - 1; j < pars.nx + BC; j++) {
        left = &(grid[i][j]);
        right = &(grid[i][j + 1]);
        solver_compute_cell_pair_flux(left, right, dt, dimension);
      }
    }
  }

#endif /* ndim */
}

void solver_prepare_flux_computation(cell *left, cell *right, int dim) {
  /* ---------------------------------------------------------------------------
   * For the WAF method, we need to pre-compute the jumps of density over
   * each wave emerging from the solution of the Riemann problem first in order
   * to be able to do flux limiting properly. So before we can compute the WAF
   * fluxes, we need to solve all Riemann problems, get all wave speeds, and
   * all densities at each cell boundary.
   *
   * cell* left:  pointer to left cell
   * cell* right: pointer to right cell
   * int dim:     along which dimension we are working
   * ---------------------------------------------------------------------------
   */

  /* first we need to solve the Riemann problem for every cell, and to store the
   * wave speeds, fluxes, and density differences of every state of the Riemann
   * solution */

  for (int k = 0; k < 4; k++) {
    gas_init_cstate(&left->riemann_fluxes[k]);
  }
  for (int k = 0; k < 3; k++) {
    left->Sk[k] = 0.;
    left->delta_q[k] = 0;
  }

  if (riemann_has_vacuum(&left->prim, &right->prim, dim)) {
    riemann_get_full_vacuum_solution_for_WAF(&left->prim, &right->prim,
                                             left->Sk, left->riemann_fluxes,
                                             left->delta_q, dim);
  } else {
#if RIEMANN == HLLC
    riemann_get_hllc_full_solution_for_WAF(&left->prim, &right->prim, left->Sk,
                                           left->riemann_fluxes, left->delta_q,
                                           dim);
#else
    riemann_get_full_solution_for_WAF(&left->prim, &right->prim, left->Sk,
                                      left->riemann_fluxes, left->delta_q, dim);
#endif
  }
}


void solver_compute_cell_pair_flux(cell *left, cell *right, float *dt,
                                   int dim) {
  /* --------------------------------------------------------------------
   * Compute the WAF flux F_{i+1/2} and store it in the left cell.
   * right: pointer to cell which stores the right state
   * dt:    current time step
   * dim:   integer along which dimension to advect. 0: x. 1: y.
   * -------------------------------------------------------------------- */

  /* get Courant numbers c_k  */
  /* ------------------------ */
  float dtdx = (*dt) / pars.dx;
  float ck[3];

  for (int k = 0; k < 3; k++) {
    ck[k] = left->Sk[k] * dtdx;
  }

  /* Find the flux limiter function phi */
  float phi[3] = {1., 1.,
                  1.}; /* set the default value for the no limiter situation */

#if LIMITER != NONE

  int i, j;
  cell *upwind = NULL;
  cell_get_ij(left, &i, &j);

  for (int k = 0; k < 3; k++) {
    /* first find the upwind cell to be able to compute r */
#if NDIM == 1
    if (ck[k] > 0.) {
      upwind = &grid[i - 1];
    } else {
      upwind = &grid[i + 1];
    }
#elif NDIM == 2
    if (ck[k] >= 0.) {
      if (dim == 0) {
        upwind = &grid[i - 1][j];
      } else if (dim == 1) {
        upwind = &grid[i][j - 1];
      }
    } else {
      if (dim == 0) {
        upwind = &grid[i + 1][j];
      } else if (dim == 1) {
        upwind = &grid[i][j + 1];
      }
    }
#endif

    /* now compute r. Check for stability when dividing first though! */
    float oneoverdq;

    if (fabsf(left->delta_q[k]) < SMALLRHO) {
      oneoverdq = 1.e6;
      if (left->delta_q[k] < 0.) {
        oneoverdq = -oneoverdq;
      }
    } else {
      oneoverdq = 1. / left->delta_q[k];
    }

    float r = upwind->delta_q[k] * oneoverdq;

    phi[k] = limiter_phi_of_r(r);
  }
#endif /* if limiter != NONE */

  /* Now do the c_k (F^(k+1) - F^(k)) sum */
  /* ------------------------------------ */

  cstate fluxsum;
  gas_init_cstate(&fluxsum);

  for (int k = 0; k < 3; k++) {

    float abscfl = fabsf(ck[k]);
    float psi = 1. - (1. - abscfl) * phi[k];
    float s = 1.;
    if (ck[k] < 0.) {
      s = -1.;
    }

    fluxsum.rho +=
        s * psi *
        (left->riemann_fluxes[k + 1].rho - left->riemann_fluxes[k].rho);
    fluxsum.rhou[0] +=
        s * psi *
        (left->riemann_fluxes[k + 1].rhou[0] - left->riemann_fluxes[k].rhou[0]);
    fluxsum.rhou[1] +=
        s * psi *
        (left->riemann_fluxes[k + 1].rhou[1] - left->riemann_fluxes[k].rhou[1]);
    fluxsum.E +=
        s * psi * (left->riemann_fluxes[k + 1].E - left->riemann_fluxes[k].E);
  }

  /* finally, compute F_i+1/2 and store it at cell i */
  /* ----------------------------------------------- */

  left->cflux.rho = 0.5 * (left->riemann_fluxes[0].rho +
                           left->riemann_fluxes[3].rho - fluxsum.rho);
  left->cflux.rhou[0] =
      0.5 * (left->riemann_fluxes[0].rhou[0] + left->riemann_fluxes[3].rhou[0] -
             fluxsum.rhou[0]);
  left->cflux.rhou[1] =
      0.5 * (left->riemann_fluxes[0].rhou[1] + left->riemann_fluxes[3].rhou[1] -
             fluxsum.rhou[1]);
  left->cflux.E =
      0.5 * (left->riemann_fluxes[0].E + left->riemann_fluxes[3].E - fluxsum.E);
}
