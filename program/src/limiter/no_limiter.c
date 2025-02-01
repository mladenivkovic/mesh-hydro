/* no limiter  */

/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */

#include "limiter.h"
#include "params.h"

extern params pars;

float limiter_phi_of_r(float r) {
  /*-------------------------------------------
   * compute the actual phi(r)
   * ------------------------------------------*/

  float phi_no_limiter = 0.5 * (1. + r); /* centered slope     */
  /* float phi_no_limiter = r;               [> upwind slope    <] */
  /* float phi_no_limiter = 1.;              [> downwind slope  <] */
  /* float phi_no_limiter = 0.;              [> no slope        <] */

  return (phi_no_limiter);
}

float limiter_xi_of_r(float r) {
  /* -----------------------------------------
   * compute the actual xi(r)
   * ----------------------------------------- */

  return (1.);
}
