/* VAN_LEER limiter  */

/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */

#include <math.h>

#include "limiter.h"
#include "params.h"

extern params pars;

/**
 * compute the actual flux limiter phi(r)
 * for the van Leer limiter
 */
float limiter_phi_of_r(float r) {

  float phi = (r + fabsf(r)) / (1. + fabsf(r));
  return (phi);
}


/**
 * compute the actual slope limiter phi(r)
 * for the van Leer limiter
 */
float limiter_xi_of_r(float r) {

  float xi = 0.;

  if (r > 0.) {
    xi        = (2. * r) / (1. + r);
    float d   = 1. - OMEGA + (1. + OMEGA) * r;
    float xiR = 2. / d;
    if (xiR < xi) { xi = xiR; }
  }

  return (xi);
}
