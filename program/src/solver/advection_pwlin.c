/* Piecewise linear advection scheme */

/* Written by Mladen Ivkovic, JAN 2020
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
  solver_get_dt(dt);
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






void solver_get_dt(float* dt){
  /* ---------------------------------------------- 
   * Computes the maximal allowable time step size
   * find max velocity present, then apply Ccfl
   * ---------------------------------------------- */

  debugmessage("Called solver_get_dt", *dt);


#if NDIM == 1
  float umax = 0;

#ifndef ADVECTION_KEEP_VELOCITY_CONSTANT
  for (int i = BC; i <= pars.nx + BC; i++){
    float uxabs = fabs(grid[i].prim.u[0]);
    if (uxabs > umax){
      umax = uxabs;
    }
  }
#else
  /* in this case, all velocities are the same and constant. */
  umax = fabs(grid[BC].prim.u[0]);
#endif

  *dt = pars.ccfl * pars.dx / umax;


#elif NDIM == 2

  float uxmax = 0;
  float uymax = 0;

#ifndef ADVECTION_KEEP_VELOCITY_CONSTANT
  for (int i = BC; i <= pars.nx + BC; i++){
    for (int j = BC; j <= pars.nx + BC; j++){
      float uxabs = fabs(grid[i][j].prim.u[0]);
      if (uxabs > uxmax){
        uxmax = uxabs;
      }
      float uyabs = fabs(grid[i][j].prim.u[1]);
      if (uyabs > uymax){
        uymax = uyabs;
      }
    }
  }
#else
  /* in this case, all velocities are the same and constant. */
  uxmax = fabs(grid[BC][BC].prim.u[0]);
  uymax = fabs(grid[BC][BC].prim.u[1]);
#endif

  float uxdx = uxmax / pars.dx; /* ux_max / dx */
  float uydy = uymax / pars.dx; /* uy_max / dy */
  
  *dt = pars.ccfl / ( uxdx + uydy );

#endif /* NDIM == 2*/

  if (pars.force_dt > 0){
    if (*dt > pars.force_dt){
      *dt = pars.force_dt;
    } else{
      throw_error("I require a smaller timestep dt=%g than force_dt=%g is demanding.",
        *dt, pars.force_dt);
    }
  }

  if (*dt <= DT_MIN) throw_error("Got weird time step? dt=%12.4e", *dt);
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
  cell *uw; /* the cell upwind of the flux at the interface of *c */
  cell *dw; /* the cell downwind of the flux at the interface of *c */

#if NDIM == 1

  for (int i = BC; i < pars.nx + BC; i++){
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

  if (dimension == 0){
    for (int i = BC; i < pars.nx + BC; i++){
      for (int j = BC; j < pars.nx + BC; j++){
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
  } else if (dimension == 1){
    for (int i = BC; i < pars.nx + BC; i++){
      for (int j = BC; j < pars.nx + BC; j++){
        c = &(grid[i][j]);
        if (c->prim.u[1] > 0) {
          dw = c;
          uw = &(grid[i][j-1]);
        } else {
          dw = &(grid[i][j+1]);
          uw = c;
        }
        solver_compute_cell_pair_flux(c, uw, dw, dt, dimension);
      }
    }

  }
#endif /* ndim */
}





void solver_compute_cell_pair_flux(cell* c, cell* uw, cell* dw, float* dt, int dim){
  /* --------------------------------------------------------------------
   * Compute the net flux for a given cell w.r.t. a specific cell pair
   * c:   pointer to cell to work with
   * uw:  the cell of the pair which is upwind
   * dw:  the cell of the cell pair which is downwind
   *
   * dim: integer along which dimension to advect. 0: x. 1: y.
   * -------------------------------------------------------------------- */

  pstate sl, sr;            /* slopes left and right of the cell */
  pstate su;                /* slope of upwind cell */
  pstate sd;                /* slope of downwind cell */

  gas_init_pstate(&sl);
  gas_init_pstate(&sr);

  limiter_get_slope_left(c, &sl, dim);
  limiter_get_slope_right(c, &sr, dim);

  /* assign the slopes to the upwind/downwind cells */
  if (uw->id < dw->id){
    su = sl;
    sd = sr;
  } else {
    su = sr;
    sd = sl;
  }

  /* Now compute the net fluxes */

  float dsu, dsd; /* dx - |v|*dt */

  if (uw->prim.u[dim] >= 0) {
    dsu = pars.dx - uw->prim.u[dim] * (*dt);
  } else {
    dsu = - pars.dx - uw->prim.u[dim] * (*dt);
  }
  if (dw->prim.u[dim] >= 0) {
    dsd = pars.dx - dw->prim.u[dim] * (*dt);
  } else {
    dsd = -pars.dx - dw->prim.u[dim] * (*dt);
  }

  c->pflux.rho += uw->prim.u[dim] * ( uw->prim.rho +  0.5 * su.rho * dsu ) -
                 dw->prim.u[dim] * ( dw->prim.rho +  0.5 * sd.rho * dsd );
#ifndef ADVECTION_KEEP_VELOCITY_CONSTANT
  c->pflux.u[0] += uw->prim.u[dim] * ( uw->prim.u[0] +  0.5 * su.u[0] * dsu ) -
                  dw->prim.u[dim] * ( dw->prim.u[0] +  0.5 * sd.u[0] * dsd );
  c->pflux.u[1] += uw->prim.u[dim] * ( uw->prim.u[1] +  0.5 * su.u[1] * dsu ) -
                  dw->prim.u[dim] * ( dw->prim.u[1] +  0.5 * sd.u[1] * dsd );
#endif
  c->pflux.p += uw->prim.u[dim] * ( uw->prim.p +  0.5 * su.p * dsu ) -
               dw->prim.u[dim] * ( dw->prim.p +  0.5 * sd.p * dsd );
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






void solver_advection_check_global_velocity(){
  /* -----------------------------------------------------
   * Check whether the velocities in the IC are constant
   * -----------------------------------------------------*/

#if NDIM == 1
  float ux = grid[BC].prim.u[0];
  for (int i = BC; i < pars.nx + BC; i++ ){
    if (grid[i].prim.u[0] != ux) {
      throw_error("The velocities are not identical everywhere. u[%d] = %12.6f; u[%d] = %12.6f\n", BC, ux, i, grid[i].prim.u[0]);
    }
  }
#elif NDIM == 2
  float ux = grid[BC][BC].prim.u[0];
  float uy = grid[BC][BC].prim.u[1];
  for (int i = BC; i < pars.nx + BC; i++ ){
    for (int j = BC; j < pars.nx + BC; j++ ){
      if (grid[i][j].prim.u[0] != ux) {
        throw_error("The velocities are not identical everywhere. ux[%d] = %12.6f; ux[%d] = %12.6f\n", BC, ux, i, grid[i][j].prim.u[0]);
      }
      if (grid[i][j].prim.u[1] != uy) {
        throw_error("The velocities are not identical everywhere. uy[%d] = %12.6f; uy[%d] = %12.6f\n", BC, uy, i, grid[i][j].prim.u[1]);
      }
    }
  }
#endif
}
