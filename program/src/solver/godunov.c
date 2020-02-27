/* Upwind Godunov Scheme */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/

#include <math.h>
#include <stdio.h>

#include "cell.h"
#include "defines.h"
#include "solver.h"
#include "params.h"
#include "utils.h"


#if NDIM == 1
extern cell *grid;
#elif NDIM == 2
extern cell **grid;
#endif

extern params pars;



void solver_init_step(){
  /* --------------------------------------------- 
   * Do everything that needs to be done before
   * we can compute the fluxes, the timestep, and
   * finally advance the simulation
   * --------------------------------------------- */
  printf("Called solver_init_step\n", *dt);
}





void solver_get_dt(float* dt){
  /* ---------------------------------------------- 
   * Computes the maximal allowable time step size
   * ---------------------------------------------- */

  printf("Called solver_get_dt with dt=%f\n", *dt);

}



void solver_compute_fluxes(){
  /* ---------------------------------------------
   * Computes the actual fluxes between cells
   * --------------------------------------------- */
   printf("Called solver_compute_fuxes\n");
}




void solver_advance_step(float* dt){
  /* ---------------------------------------------
   * Integrate the equations for one time step
   * --------------------------------------------- */

  printf("Called solver_advance_step with dt=%f\n", *dt);
}
