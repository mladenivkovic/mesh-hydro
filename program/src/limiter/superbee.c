/* SUPERBEE limiter  */

/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */

#include "cell.h"
#include "defines.h"
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
  /* ---------------------------------------
   * compute the actual flux limiter phi(r)
   * for the superbee limiter
   * --------------------------------------- */

  float max = 0.0;

  /* min(1, 2r) */
  float temp = 1.0;
  if (temp > 2 * r)
    temp = 2 * r;
  if (temp > max)
    max = temp;

  /* min(2, r) */
  temp = 2.0;
  if (temp > r)
    temp = r;
  if (temp > max)
    max = temp;
  return (max);
}

float limiter_xi_of_r(float r) {
  /* ---------------------------------------
   * compute the actual slope limiter xi(r)
   * for the superbee limiter
   * --------------------------------------- */
  float xi = 0;
  if (r > 0.)
    xi = 2. * r;
  if (r > 0.5)
    xi = 1.0;
  if (r > 1.) {
    float d = 1. - OMEGA + (1. + OMEGA) * r;
    float xiR = 2. / d;

    /* xi = min(xiR, r) */
    xi = r;
    if (xiR < xi)
      xi = xiR;

    /* xi = min(xi, 2) */
    if (xi > 2.)
      xi = 2.;
  }
  return (xi);
}
