/* VAN_LEER limiter  */


/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */

#include "gas.h"
#include "cell.h"
#include "limiter.h"
#include "params.h"

#include <math.h>
#include <stdio.h>


#if NDIM == 1
extern cell *grid;
#elif NDIM == 2
extern cell **grid;
#endif


extern params pars;






float limiter_phi_of_r(float r){
  /* -----------------------------
   * compute the actual phi(r) 
   * ----------------------------- */

  return(limiter_vanleer(r));
}





float limiter_vanleer(float r){
  /* -----------------------------------------------------
   * Computes the phi(r) for the van Leer slope limiter
   * ----------------------------------------------------- */

  float phi = (r + fabs(r))/(1 + fabs(r));
  return(phi);
}
