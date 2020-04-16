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
  /* -----------------------------------------
   * compute the actual flux limiter phi(r) 
   * for the van Leer limiter
   * ----------------------------------------- */

  float phi = (r + fabs(r))/(1 + fabs(r));
  return(phi);
}





float limiter_xi_of_r(float r){
  /* -----------------------------------------
   * compute the actual slope limiter phi(r) 
   * for the van Leer limiter
   * ----------------------------------------- */

  float xi = 0;

  if (r > 0.) {
    xi = (2. * r)/(1. + r);
    float d = 1. - OMEGA + (1. + OMEGA)*r;
    float xiR = 2./d;
    if (xiR < xi) xi = xiR;
  }

  return(xi);
}
