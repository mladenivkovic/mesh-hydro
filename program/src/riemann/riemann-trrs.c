/* Two Rarefaction Approximate Riemann Solver */

/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */

#include "gas.h"
#include "params.h"
#include "riemann.h"
#include "utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern params pars;

void riemann_compute_star_states(pstate *left, pstate *right, float *pstar,
                                 float *ustar, int dim) {
  /* ------------------------------------------------------------------------------------------
   * computes the star region pressure and velocity given the left and right
   *pstates.
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * float* pstar:    where pressure of star region will be written
   * float* ustar:    where velocity of star region will be written
   * int dim:         which fluid velocity dimensionto use. 0: x, 1: y
   *------------------------------------------------------------------------------------------*/

  float aL = gas_soundspeed(left);
  float aLinv = 1. / aL;
  float aR = gas_soundspeed(right);
  float aRinv = 1. / aR;

  float pLRbeta = pow(left->p / right->p, BETA);

  *ustar = ((pLRbeta - 1.) / GM1HALF + left->u[dim] * aLinv * pLRbeta +
            right->u[dim] * aRinv) /
           (aRinv + aLinv * pLRbeta);

  *pstar =
      0.5 * (right->p * pow((1. + aRinv * GM1HALF * (*ustar - right->u[dim])),
                            1. / BETA) +
             left->p * pow((1. + aLinv * GM1HALF * (left->u[dim] - *ustar)),
                           1. / BETA));

  if (*pstar < SMALLP)
    *pstar = SMALLP;

  debugmessage("Got pstar = %12.8f, ustar = %12.8f", *pstar, *ustar);
}
