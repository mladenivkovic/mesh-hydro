/* Weighted Average Flux advection scheme 

 * Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/

#include <math.h>
#include <stdio.h>

#include "cell.h"
#include "defines.h"
#include "solver.h"
#include "io.h"
#include "limiter.h"
#include "params.h"
#include "utils.h"


#if NDIM == 1
extern cell *grid;
#elif NDIM == 2
extern cell **grid;
#endif

extern params pars;





void solver_step(float *t, float* dt, int step, int* write_output){
  /* -------------------------------------------------------
   * Main routine for the actual hydro step
   * ------------------------------------------------------- */

  solver_init_step();
  solver_get_advection_dt(dt);
  /* check this here in case you need to limit time step for output */
  *write_output = io_is_output_step(*t, dt, step); 

#if NDIM == 1

  solver_compute_fluxes(dt, /*dimension =*/0);
  solver_advance_step(dt);

#elif NDIM == 2

  int dimension = step % 2; /* gives 0 or 1, switching each step */
  solver_compute_fluxes(dt, dimension);
  solver_advance_step(dt);

  dimension = (dimension + 1) % 2; /* 1 -> 0 or 0 -> 1 */
  solver_init_step();
  solver_compute_fluxes(dt, dimension);
  solver_advance_step(dt);

#endif
}




void solver_init_step(){
  /* --------------------------------------------- 
   * Do everything that needs to be done before
   * we can compute the fluxes, the timestep, and
   * finally advance the simulation
   * --------------------------------------------- */

  debugmessage("Called solver_init_step");
  cell_reset_fluxes();
  cell_set_boundary();
}






void solver_compute_fluxes(float* dt, int dimension){
  /* ------------------------------------------------------
   * Computes the actual *net* fluxes between cells
   * Here we compute F^{n+1/2}_{i-1/2} - F^{n+1/2}_{i+1/2}
   * and store it in cell.pflux
   * dimension: 0 for x, 1 for y. Only used with Strang
   * splitting.
   * ------------------------------------------------------ */

  debugmessage("Called solver_compute_fluxes; dimension = %d", dimension);

  cell *c; /* this cell */
  cell *n; /* the neighbour we're dealing with, i.e. i+1 or j+1 */

#if NDIM == 1

  /* computes F_{i+1/2}, so we need to start at index BC-1:
   * F_{i+1/2} = F_{(i+1) - 1/2}*/
  for (int i = BC-1; i < pars.nx + BC; i++){
    c = &grid[i];
    n = &grid[i+1];
    solver_compute_cell_pair_flux(c, n, dt, /*dimension=*/0);
  }

#elif NDIM == 2

  if (dimension == 0){
    /* computes F_{i+1/2}, so we need to start at index BC-1:
     * F_{i+1/2} = F_{(i+1) - 1/2}*/
    for (int i = BC-1; i < pars.nx + BC; i++){
      for (int j = BC-1; j < pars.nx + BC; j++){
        c = &(grid[i][j]);
        n = &(grid[i+1][j]);
        solver_compute_cell_pair_flux(c, n, dt, dimension);
      }
    }
  } else if (dimension == 1){
    /* computes F_{i+1/2}, so we need to start at index BC-1:
     * F_{i+1/2} = F_{(i+1) - 1/2}*/
    for (int i = BC-1; i < pars.nx + BC; i++){
      for (int j = BC-1; j < pars.nx + BC; j++){
        c = &(grid[i][j]);
        n = &(grid[i][j+1]);
        solver_compute_cell_pair_flux(c, n, dt, dimension);
      }
    }

  }
#endif /* ndim */
}





void solver_compute_cell_pair_flux(cell* c, cell* n, float* dt, int dim){
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
 
  float cfl = (*dt) / pars.dx * c->prim.u[dim];
  float flux = 0.;

  pstate psi;
  gas_init_pstate(&psi);
  limiter_get_psi(c, &psi, cfl, dim);

  pstate phi;
  gas_init_pstate(&phi);
  phi.rho = 1. - (1. - fabs(cfl)) * psi.rho;
  phi.u[0] = 1. - (1. - fabs(cfl)) * psi.u[0];
  phi.u[1] = 1. - (1. - fabs(cfl)) * psi.u[1];
  phi.p = 1. - (1. - fabs(cfl)) * psi.p;



  float vel = c->prim.u[dim];
  float s = 1.; /* sign(velocity) */
  if (vel <= 0) s = -1.;


  flux = 0.5 * (1. + s*phi.rho) * vel * c->prim.rho +
         0.5 * (1. - s*phi.rho) * vel * n->prim.rho;
  c->pflux.rho -= flux;
  n->pflux.rho += flux;
#ifndef ADVECTION_KEEP_VELOCITY_CONSTANT
  flux = 0.5 * (1. + s*phi.u[0]) * vel * c->prim.u[0] +
         0.5 * (1. - s*phi.u[0]) * vel * * n->prim.u[0];
  c->pflux.u[0] -= flux;
  n->pflux.u[0] += flux;
  flux = 0.5 * (1. + s*phi.u[1]) * vel * c->prim.u[1] +
         0.5 * (1. - s*phi.u[1]) * vel * n->prim.u[1];
  c->pflux.u[1] -= flux;
  n->pflux.u[1] += flux;
#endif
  flux = 0.5 * (1. + s*phi.p) * vel * c->prim.p +
         0.5 * (1. - s*phi.p) * vel * n->prim.p;
  c->pflux.p -= flux;
  n->pflux.p += flux;
}





void solver_advance_step(float* dt){
  /* ---------------------------------------------
   * Integrate the equations for one time step
   * --------------------------------------------- */

  debugmessage("Called solver_advance_step with dt = %f", *dt);

  float dtdx = *dt / pars.dx;

#if NDIM == 1
  for (int i = BC; i < pars.nx + BC; i++){
    solver_update_state(&(grid[i]), dtdx);
  }
#elif NDIM == 2
  for (int i = BC; i < pars.nx + BC; i++){
    for (int j = BC; j < pars.nx + BC; j++){
      solver_update_state(&(grid[i][j]), dtdx);
    }
  }
#endif
}





void solver_update_state(cell *c, float dtdx){
  /* ------------------------------------------------------
   * Update the state using the fluxes in the cell and dt
   * dtdx: dt / dx
   * ------------------------------------------------------ */

  c->prim.rho = c->prim.rho + dtdx * c->pflux.rho;
#ifndef ADVECTION_KEEP_VELOCITY_CONSTANT
  c->prim.u[0] = c->prim.u[0] + dtdx * c->pflux.u[0];
#if NDIM >= 2
  c->prim.u[1] = c->prim.u[1] + dtdx * c->pflux.u[1];
#endif
#endif
  c->prim.p = c->prim.p + dtdx * c->pflux.p;
}
