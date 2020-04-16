/* monotonized central difference limiter  */


/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */

#include "gas.h"
#include "cell.h"
#include "limiter.h"
#include "params.h"
#include "utils.h"

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

  return(limiter_mc(r));
}



float limiter_xi_of_r(float r){
  /* ----------------------------
   * compute slope limiter xi(r)
   * NOT IMPLEMENTED FOR MC!!!
   * ---------------------------- */
  throw_error("The slope limiter \\xi(r) is not implemented for the MC limiter. \nPlease try something else");
  return(1.);
}



float limiter_mc(float r){
  /* -----------------------------------------------------
   * Computes the phi(r) for the MC slope limiter
   * ----------------------------------------------------- */

  float min = 0.0;
  float max = 0.0;

  /* min((1+r)/2, 2) */
  float temp = 0.5*(1+r);
  if (temp > 2.0) temp = 2.0;
  min = temp;

  /* min( min((1+r)/2, 2), 2r) */
  temp = 2 * r;
  if (temp > min) temp = min;       
  min = temp;
  
  /* max( 0, min(...) ) */
  if (min > max) max = min;

  return(max);
}
