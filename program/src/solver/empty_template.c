/* WRITE SOMETHING HERE */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/


#include <math.h>
#include <stdio.h>

#include "cell.h"
#include "defines.h"
#include "solver.h"
#include "io.h"
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
    solver_compute_fluxes();

    /* check this here in case you need to limit time step for output */
    *write_output = io_is_output_step(*t, dt, step); 

    solver_advance_step(dt);
}



void solver_init_step(){
  /* --------------------------------------------- 
   * Do everything that needs to be done before
   * we can compute the fluxes, the timestep, and
   * finally advance the simulation
   * --------------------------------------------- */
  debugmessage("Called solver_init_step");
}





void solver_get_dt(MYFLOAT* dt){
  /* ---------------------------------------------- 
   * Computes the maximal allowable time step size
   * ---------------------------------------------- */

  debugmessage("Called solver_get_dt with dt=%f", *dt);
  *dt = 0.12345;

}



void solver_compute_fluxes(){
  /* ---------------------------------------------
   * Computes the actual fluxes between cells
   * --------------------------------------------- */
   debugmessage("Called solver_compute_fluxes");
}




void solver_advance_step(MYFLOAT* dt){
  /* ---------------------------------------------
   * Integrate the equations for one time step
   * --------------------------------------------- */

  debugmessage("Called solver_advance_step with dt=%f", *dt);
}
