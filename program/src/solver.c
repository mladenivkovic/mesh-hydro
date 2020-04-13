/* some method independent routines *
 * Written by Mladen Ivkovic, MAR 2020
 * mladen.ivkovic@hotmail.com           */

#include <math.h>

#include "cell.h"
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






void solver_get_advection_dt(float* dt){
  /* ---------------------------------------------- 
   * Computes the maximal allowable time step size
   * find max velocity present, then apply 
   * the Courant number.
   *
   * Intended for advection methods.
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







void solver_get_hydro_dt(float* dt, int step){
  /* ---------------------------------------------- 
   * Computes the maximal allowable time step size
   * find max velocity present, then apply Courant
   * number to reduce it.
   * ---------------------------------------------- */

  debugmessage("Called solver_get_hydro_dt", *dt);


#if NDIM == 1
  float umax = 0.;

  for (int i = BC; i < pars.nx + BC; i++){
    float uxabs = fabs(grid[i].prim.u[0]);
    float a = gas_soundspeed(&grid[i].prim);
    float S = uxabs + a;
    if (S > umax){ umax = S; }
  }

  *dt = pars.ccfl * pars.dx / umax;

#elif NDIM == 2

  float uxmax = 0;
  float uymax = 0;

  for (int i = BC; i < pars.nx + BC; i++){
    for (int j = BC; j < pars.nx + BC; j++){
      float uxabs = fabs(grid[i][j].prim.u[0]);
      float a = gas_soundspeed(&grid[i][j].prim);
      float S = uxabs + a;
      if (S > uxmax){ uxmax = S; }
      float uyabs = fabs(grid[i][j].prim.u[1]);
      S = uyabs + a;
      if (S > uymax){ uymax = S; }
    }
  }

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

  /* sometimes there might be trouble with sharp discontinuities at the beginning, 
   * so reduce the timestep for the first few steps */
  if (step <=5) *dt *= 0.2;

  /* safety check */
  if (*dt <= DT_MIN) throw_error("Got weird time step? dt=%12.4e", *dt);
}


