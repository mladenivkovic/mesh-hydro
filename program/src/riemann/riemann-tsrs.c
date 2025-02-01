/* Two Shock Riemann Solver */

/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */

#include "gas.h"
#include "params.h"
#include "riemann.h"
/* #include "utils.h" */

#include <math.h>

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
   * int dimension:   which fluid velocity dimensionto use. 0: x, 1: y
   *------------------------------------------------------------------------------------------*/

  float AL = 2. / GP1 / left->rho;
  float AR = 2. / GP1 / right->rho;
  float BL = GM1OGP1 * left->p;
  float BR = GM1OGP1 * right->p;
  float aL = gas_soundspeed(left);
  float aR = gas_soundspeed(right);

  float delta_u = right->u[dim] - left->u[dim];

  /* Find initial guess for star pressure */
  float pguess = 0.5 * (left->p + right->p) -
                 0.125 * delta_u * (left->rho + right->rho) * (aL + aR);
  if (pguess < EPSILON_ITER) {
    pguess = EPSILON_ITER;
  }

  float gL = sqrtf(AL / (pguess + BL));
  float gR = sqrtf(AR / (pguess + BR));
  *pstar = (gL * left->p + gR * right->p - delta_u) / (gL + gR);
  if (SMALLP > *pstar) {
    *pstar = SMALLP;
  }

  *ustar = 0.5 * (right->u[dim] + left->u[dim] + (*pstar - right->p) * gR -
                  (*pstar - left->p) * gL);

  /* debugmessage("Got pstar = %12.8f, ustar = %12.8f", *pstar, *ustar); */
}
