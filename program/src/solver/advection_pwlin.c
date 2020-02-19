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





void solver_step(MYFLOAT *t, MYFLOAT* dt, int step, int* write_output){
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
#ifndef STRANG_SPLITTING
  solver_compute_fluxes(/*dimension =*/0);
  solver_advance_step(dt);
#else
  int dimension = step % 2; /* gives 0 or 1, switching each step */
  solver_compute_fluxes(dt, dimension);
  solver_advance_step(dt);
  cell_reset_fluxes();
  dimension = (dimension + 1) % 2; /* 1 -> 0 or 0 -> 1 */
  solver_compute_fluxes(dt, dimension);
  solver_advance_step(dt);
  /* cell_reset_fluxes(); */ /* will be done in solver_init_step */
#endif
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






void solver_get_dt(MYFLOAT* dt){
  /* ---------------------------------------------- 
   * Computes the maximal allowable time step size
   * find max velocity present, then apply Ccfl
   * ---------------------------------------------- */

  debugmessage("Called solver_get_dt", *dt);


#if NDIM == 1
  MYFLOAT umax = 0;

#ifndef ADVECTION_KEEP_VELOCITY_CONSTANT
  for (int i = BC; i <= pars.nx + BC; i++){
    MYFLOAT uxabs = fabs(grid[i].prim.ux);
    if (uxabs > umax){
      umax = uxabs;
    }
  }
#else
  /* in this case, all velocities are the same and constant. */
  umax = fabs(grid[BC].prim.ux);
#endif

  *dt = pars.ccfl * pars.dx / umax;


#elif NDIM == 2

  MYFLOAT uxmax = 0;
  MYFLOAT uymax = 0;

#ifndef ADVECTION_KEEP_VELOCITY_CONSTANT
  for (int i = BC; i <= pars.nx + BC; i++){
    for (int j = BC; j <= pars.nx + BC; j++){
      MYFLOAT uxabs = fabs(grid[i][j].prim.ux);
      if (uxabs > uxmax){
        uxmax = uxabs;
      }
      MYFLOAT uyabs = fabs(grid[i][j].prim.uy);
      if (uyabs > uymax){
        uymax = uyabs;
      }
    }
  }
#else
  /* in this case, all velocities are the same and constant. */
  uxmax = fabs(grid[BC][BC].prim.ux);
  uymax = fabs(grid[BC][BC].prim.uy);
#endif

  MYFLOAT uxdx = uxmax / pars.dx; /* ux_max / dx */
  MYFLOAT uydy = uymax / pars.dx; /* uy_max / dy */
  
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

  if (*dt <= 0.0) throw_error("Got weird time step? dt=%12.8f");
}






void solver_compute_fluxes(MYFLOAT* dt, int dimension){
  /* ------------------------------------------------------
   * Computes the actual *net* fluxes between cells
   * Here we compute F^{n+1/2}_{i-1/2} - F^{n+1/2}_{i+1/2}
   * and store it in cell.flux
   * dimension: 0 for x, 1 for y. Only used with Strang
   * splitting.
   * ------------------------------------------------------ */

  debugmessage("Called solver_compute_fluxes");

  cell *c; /* this cell */
  cell *uw; /* the cell upwind of the flux at the interface of *c */
  cell *dw; /* the cell downwind of the flux at the interface of *c */

#if NDIM == 1

  for (int i = BC; i < pars.nx + BC; i++){
    c = &grid[i];
    if (c->prim.ux > 0) { /* we do upwind differencing */
      dw = c;
      uw = &(grid[i - 1]);
    } else {
      dw = &(grid[i + 1]);
      uw = c;
    }
    solver_compute_cell_pair_flux(c, uw, dw, dt, /*dimension=*/0);
  }

#elif NDIM == 2
#ifdef STRANG_SPLITTING

  if (dimension == 0){
    for (int i = BC; i < pars.nx + BC; i++){
      for (int j = BC; j < pars.nx + BC; j++){
        c = &(grid[i][j]);
        if (c->prim.ux > 0) {
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
        if (c->prim.uy > 0) {
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
#else
  for (int i = BC; i < pars.nx + BC; i++){
    for (int j = BC; j < pars.nx + BC; j++){
      c = &(grid[i][j]);
      if (c->prim.ux > 0) {
        dw = c;
        uw = &(grid[i - 1][j]);
      } else {
        dw = &(grid[i + 1][j]);
        uw = c;
      }
      solver_compute_cell_pair_flux(c, uw, dw, /*dimension=*/0);
      if (c->prim.uy > 0) {
        dw = c;
        uw = &(grid[i][j-1]);
      } else {
        dw = &(grid[i][j+1]);
        uw = c;
      }
      solver_compute_cell_pair_flux(c, uw, dw, /*dimension=*/1);
    }
  }
#endif /* strang splitting */
#endif /* ndim */
}





void solver_compute_cell_pair_flux(cell* c, cell* uw, cell* dw, MYFLOAT* dt, int dimension){
  /* --------------------------------------------------------------------
   * Compute the net flux for a given cell w.r.t. a specific cell pair
   * c:   pointer to cell to work with
   * uw:  the cell of the pair which is upwind
   * dw:  the cell of the cell pair which is downwind
   *
   * dimension: integer along which dimension to advect. 0: x. 1: y.
   * -------------------------------------------------------------------- */

  pstate sl = {0, 0, 0, 0};
  pstate sr = {0, 0, 0, 0}; /* slopes left and right of the cell */
  pstate su;                /* slope of upwind cell */
  pstate sd;                /* slope of downwind cell */

  limiter_get_slope_left(c, &sl, dimension);
  limiter_get_slope_right(c, &sr, dimension);

  /* assign the slopes to the upwind/downwind cells */
  if (uw->id < dw->id){
    su = sl;
    sd = sr;
  } else {
    su = sr;
    sd = sl;
  }

  if (dimension == 0){
  
    /* Now compute the net fluxes */

    MYFLOAT dsu, dsd; /* dx - |v|*dt */

    if (uw->prim.ux >= 0) {
      dsu = pars.dx - uw->prim.ux * (*dt);
    } else {
      dsu = - pars.dx - uw->prim.ux * (*dt);
    }
    if (dw->prim.ux >= 0) {
      dsd = pars.dx - dw->prim.ux * (*dt);
    } else {
      dsd = -pars.dx - dw->prim.ux * (*dt);
    }

    c->flux.rho += uw->prim.ux *( uw->prim.rho +  0.5 * su.rho * dsu ) -
                   dw->prim.ux *( dw->prim.rho +  0.5 * sd.rho * dsd );
#ifndef ADVECTION_KEEP_VELOCITY_CONSTANT
    c->flux.ux += uw->prim.ux *( uw->prim.ux +  0.5 * su.ux * dsu ) -
                   dw->prim.ux *( dw->prim.ux +  0.5 * sd.ux * dsd );
    c->flux.uy += uw->prim.ux *( uw->prim.uy +  0.5 * su.uy * dsu ) -
                   dw->prim.ux *( dw->prim.uy +  0.5 * sd.uy * dsd );
#endif
    c->flux.p += uw->prim.ux *( uw->prim.p +  0.5 * su.p * dsu ) -
                   dw->prim.ux *( dw->prim.p +  0.5 * sd.p * dsd );


  } else if (dimension == 1){

    /* Now compute the net fluxes */

    MYFLOAT dsu, dsd; /* dx - |v|*dt */

    if (uw->prim.uy >= 0) {
      dsu = pars.dx - uw->prim.uy * (*dt);
    } else {
      dsu = - pars.dx - uw->prim.uy * (*dt);
    }
    if (dw->prim.uy >= 0) {
      dsd = pars.dx - dw->prim.uy * (*dt);
    } else {
      dsd = -pars.dx - dw->prim.uy * (*dt);
    }

    c->flux.rho += uw->prim.uy *( uw->prim.rho +  0.5 * su.rho * dsu ) -
                   dw->prim.uy *( dw->prim.rho +  0.5 * sd.rho * dsd );
#ifndef ADVECTION_KEEP_VELOCITY_CONSTANT
    c->flux.ux += uw->prim.uy *( uw->prim.ux +  0.5 * su.ux * dsu ) -
                   dw->prim.uy *( dw->prim.ux +  0.5 * sd.ux * dsd );
    c->flux.uy += uw->prim.uy *( uw->prim.uy +  0.5 * su.uy * dsu ) -
                   dw->prim.uy *( dw->prim.uy +  0.5 * sd.uy * dsd );
#endif
    c->flux.p += uw->prim.uy *( uw->prim.p +  0.5 * su.p * dsu ) -
                   dw->prim.uy *( dw->prim.p +  0.5 * sd.p * dsd );

  }
}








void solver_advance_step(MYFLOAT* dt){
  /* ---------------------------------------------
   * Integrate the equations for one time step
   * --------------------------------------------- */

  debugmessage("Called solver_advance_step with dt = %f", *dt);

  MYFLOAT dtdx = *dt / pars.dx;

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





void solver_update_state(cell *c, MYFLOAT dtdx){
    /* ------------------------------------------------------
     * Update the state using the fluxes in the cell and dt
     * dtdx: dt / dx
     * ------------------------------------------------------ */

    c->prim.rho = c->prim.rho + dtdx * c->flux.rho;
#ifndef ADVECTION_KEEP_VELOCITY_CONSTANT
    c->prim.ux = c->prim.ux + dtdx * c->flux.ux;
#if NDIM >= 2
    c->prim.uy = c->prim.uy + dtdx * c->flux.uy;
#endif
#endif
    c->prim.p = c->prim.p + dtdx * c->flux.p;
}






void solver_advection_check_global_velocity(){
  /* -----------------------------------------------------
   * Check whether the velocities in the IC are constant
   * -----------------------------------------------------*/

#if NDIM == 1
  MYFLOAT ux = grid[BC].prim.ux;
  for (int i = BC; i < pars.nx + BC; i++ ){
    if (grid[i].prim.ux != ux) {
      throw_error("The velocities are not identical everywhere. u[%d] = %12.6f; u[%d] = %12.6f\n", BC, ux, i, grid[i].prim.ux);
    }
  }
#elif NDIM == 2
  MYFLOAT ux = grid[BC][BC].prim.ux;
  MYFLOAT uy = grid[BC][BC].prim.uy;
  for (int i = BC; i < pars.nx + BC; i++ ){
    for (int j = BC; j < pars.nx + BC; j++ ){
      if (grid[i][j].prim.ux != ux) {
        throw_error("The velocities are not identical everywhere. ux[%d] = %12.6f; ux[%d] = %12.6f\n", BC, ux, i, grid[i][j].prim.ux);
      }
      if (grid[i][j].prim.uy != uy) {
        throw_error("The velocities are not identical everywhere. uy[%d] = %12.6f; uy[%d] = %12.6f\n", BC, uy, i, grid[i][j].prim.uy);
      }
    }
  }
#endif
}
