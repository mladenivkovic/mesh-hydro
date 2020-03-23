/* SUPERBEE limiter  */


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

  return(limiter_superbee(r));
}




float limiter_superbee(float r){
  /* -----------------------------------------------------
   * Computes the phi(r) for the superbee slope limiter
   * ----------------------------------------------------- */

  float max = 0.0;

  /* min(1, 2r) */
  float temp = 1.0;
  if (temp > 2*r) temp = 2*r;
  if (temp > max) max = temp;

  /* min(2, r) */
  temp = 2.0;
  if (temp > r) temp = r;       
  if (temp > max) max = temp;
  return(max);
}
