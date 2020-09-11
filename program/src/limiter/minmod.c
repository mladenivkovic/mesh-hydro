/* MINMOD limiter  */

/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */

#include "cell.h"
#include "gas.h"
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

float limiter_phi_of_r(float r) {
  /*--------------------------------------------
   * compute the actual flux limiter phi(r) for
   * the minmod limiter
   * -------------------------------------------*/
  return (minmod(1., r));
}

float limiter_xi_of_r(float r) {
  /* -------------------------------------------
   * Compute the actual slope limiter xi(r) for
   * the minmod limiter
   * -------------------------------------------*/

  float xi = 0;
  if (r > 0.)
    xi = r;
  if (r > 1.) {
    float d = 1. - OMEGA + (1. + OMEGA) * r;
    float xiR = 2. / d;

    /* xi = min(1, xiR) */
    xi = 1.;
    if (xiR < xi)
      xi = xiR;
  }
  return (xi);
}

float minmod(float a, float b) {
  /* ----------------------------------------------
   * Computes the minmod() operator on a and b
   * ---------------------------------------------- */

  if (a * b <= 0.0)
    return (0.0);

  if (fabs(a) < fabs(b)) {
    return (a);
  } else {
    return (b);
  }
}
