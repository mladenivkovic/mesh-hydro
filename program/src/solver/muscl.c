/* Weighted Average Flux Hydro Scheme  */

/* Written by Mladen Ivkovic, APR 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/

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

/**
 * Main routine for the actual hydro step
 */
void solver_step(const float* t, float* dt, int step, int* write_output) {

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
  debugmessage(
    "Computing source term ODE for half timestep after hydro step "
    "has finished"
  );
  sources_update_state(dthalf);
#endif
}


/**
 * Do everything that needs to be done before
 * we can compute the fluxes, the timestep, and
 * finally advance the simulation
 */
void solver_init_step(void) {
  debugmessage("Called solver_init_step");
  cell_reset_fluxes();
  cell_get_pstates_from_cstates();
  cell_set_boundary();
}


/**
 * Compute the flux F_{i+1/2} (or G_{i+1/2} if dimension == 1)
 * and store it in cell i.
 * BUT: Start with the first boundary cell, such that even the
 * first real cell can access F_{i-1/2} by accessing the
 * neighbour at i-1
 * However, we first need to pre-compute the intermediate
 * boundary extrapolated states for every cell individually.
 */
void solver_compute_fluxes(float* dt, int dimension) {

  debugmessage("Called solver_compute_fluxes; dimension = %d", dimension);

  float dthalf = 0.5 * (*dt);
  cell* left;  /* this cell */
  cell* right; /* the right neighbour */

#if NDIM == 1
  /* compute intermediate boundary extrapolated states first */
  for (int i = BC - 1; i < pars.nx + BC + 1; i++) {
    solver_prepare_flux_computation(&grid[i], dthalf, /*dim=*/0);
  }

  /* now update states for this dimension */
  for (int i = BC - 1; i < pars.nx + BC; i++) {
    left  = &grid[i];
    right = &grid[i + 1];
    solver_compute_cell_pair_flux(left, right, dt, /*dim=*/0);
  }

#elif NDIM == 2

  /* compute intermediate boundary extrapolated states first */
  for (int i = BC - 1; i < pars.nx + BC + 1; i++) {
    for (int j = BC - 1; j < pars.nx + BC + 1; j++) {
      solver_prepare_flux_computation(&grid[i][j], dthalf, dimension);
    }
  }

  /* now update states for this dimension */
  if (dimension == 0) {
    for (int i = BC - 1; i < pars.nx + BC; i++) {
      for (int j = BC - 1; j < pars.nx + BC; j++) {
        left  = &(grid[i][j]);
        right = &(grid[i + 1][j]);
        solver_compute_cell_pair_flux(left, right, dt, dimension);
      }
    }
  } else if (dimension == 1) {
    for (int i = BC - 1; i < pars.nx + BC; i++) {
      for (int j = BC - 1; j < pars.nx + BC; j++) {
        left  = &(grid[i][j]);
        right = &(grid[i][j + 1]);
        solver_compute_cell_pair_flux(left, right, dt, dimension);
      }
    }
  }

#endif /* ndim */
}

/**
 * For the MUSCL-Hancock scheme, we need to first compute the slopes for each
 * conserved variable and each cell, and then compute the updated boundary
 * extrapolated values. Only then can we correctly compute the intercell
 * fluxes.
 * This function first computes the fluxes, and then computes the updated
 * intermediate state for each cell, and stores them in the cell.
 *
 * @param c:       pointer to cell to work with
 * @param dthalf:  Delta t / 2
 * @param dim:     along which dimension we are working
 */
void solver_prepare_flux_computation(cell* c, float dthalf, int dim) {

  /* first get the slope */
  cstate slope;
  gas_init_cstate(&slope);
  limiter_get_limited_slope(c, &slope, dim);

  /* compute fluxes */
  cstate UL;
  gas_init_cstate(&UL);
  UL.rho     = c->cons.rho - 0.5 * slope.rho;
  UL.rhou[0] = c->cons.rhou[0] - 0.5 * slope.rhou[0];
  UL.rhou[1] = c->cons.rhou[1] - 0.5 * slope.rhou[1];
  UL.E       = c->cons.E - 0.5 * slope.E;

  cstate FL;
  gas_init_cstate(&FL);
  gas_get_cflux_from_cstate(&UL, &FL, dim);

  cstate UR;
  gas_init_cstate(&UR);
  UR.rho     = c->cons.rho + 0.5 * slope.rho;
  UR.rhou[0] = c->cons.rhou[0] + 0.5 * slope.rhou[0];
  UR.rhou[1] = c->cons.rhou[1] + 0.5 * slope.rhou[1];
  UR.E       = c->cons.E + 0.5 * slope.E;

  cstate FR;
  gas_init_cstate(&FR);
  gas_get_cflux_from_cstate(&UR, &FR, dim);

  /* now compute intermediate cell state */
  float dtdxhalf = dthalf / pars.dx;

  c->ULmid.rho = c->cons.rho + dtdxhalf * (FL.rho - FR.rho) - 0.5 * slope.rho;
  c->ULmid.rhou[0] = c->cons.rhou[0] + dtdxhalf * (FL.rhou[0] - FR.rhou[0])
                     - 0.5 * slope.rhou[0];
  c->ULmid.rhou[1] = c->cons.rhou[1] + dtdxhalf * (FL.rhou[1] - FR.rhou[1])
                     - 0.5 * slope.rhou[1];
  c->ULmid.E = c->cons.E + dtdxhalf * (FL.E - FR.E) - 0.5 * slope.E;

  c->URmid.rho = c->cons.rho + dtdxhalf * (FL.rho - FR.rho) + 0.5 * slope.rho;
  c->URmid.rhou[0] = c->cons.rhou[0] + dtdxhalf * (FL.rhou[0] - FR.rhou[0])
                     + 0.5 * slope.rhou[0];
  c->URmid.rhou[1] = c->cons.rhou[1] + dtdxhalf * (FL.rhou[1] - FR.rhou[1])
                     + 0.5 * slope.rhou[1];
  c->URmid.E = c->cons.E + dtdxhalf * (FL.E - FR.E) + 0.5 * slope.E;
}


/**
 * Compute the flux F_{i+1/2} for a given cell w.r.t. a specific cell pair
 *
 * Here, we just solve the Riemann problem  with U_L = U^R_{i,BEXT},
 * U_R = U^L_{i+1, BEXT}, where
 *    U^R_{i,BEXT} is the intermediate right extrapolated boundary value of
 * cell i U^L_{i+1, BEXT}  is the intermediate left extrapolated boundary
 * value of cell i+1 and then sample the solution at x = x/t = 0, because that
 * is where we set the initial boundary in the local coordinate system between
 * the left and right cell.
 *
 * @param left:  pointer to cell which stores the left state
 * @param right: pointer to cell which stores the right state
 * @param dt:   current time step
 * @param dim:     integer along which dimension to advect. 0: x. 1: y.
 */
void solver_compute_cell_pair_flux(
  cell* left, cell* right, const float* dt, int dim
) {

  pstate WL;
  gas_init_pstate(&WL);
  gas_cons_to_prim(&left->URmid, &WL);

  pstate WR;
  gas_init_pstate(&WR);
  gas_cons_to_prim(&right->ULmid, &WR);

#if RIEMANN == HLLC
  /* the HLLC solver gives us the flux directly. */
  riemann_solve_hllc(&WL, &WR, &left->cflux, /*xovert=*/0.0, dim);
#else
  pstate solution;
  gas_init_pstate(&solution);

  /* solve the Riemann problem */
  riemann_solve(&WL, &WR, &solution, /*xovert=*/0.0, dim);

  /* from the primitive states, compute and store F_{i+1/2} */
  gas_get_cflux_from_pstate(&solution, &left->cflux, dim);
#endif
}
