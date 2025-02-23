/* Exact Riemann Solver */

/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */

#include "gas.h"
#include "params.h"
#include "riemann.h"
/* #include "utils.h" */

#include <math.h>
#include <stdio.h>

extern params pars;

void riemann_compute_star_states(
  pstate* left, pstate* right, float* pstar, float* ustar, int dim
) {
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

  float AL = 2. / GP1 / left->rho;
  float AR = 2. / GP1 / right->rho;
  float BL = GM1OGP1 * left->p;
  float BR = GM1OGP1 * right->p;
  float aL = gas_soundspeed(left);
  float aR = gas_soundspeed(right);

  float delta_u = right->u[dim] - left->u[dim];

  float pguess, pold;

  /* Find initial guess for star pressure */
  float ppv = 0.5 * (left->p + right->p)
              - 0.125 * delta_u * (left->rho + right->rho) * (aL + aR);
  pguess = ppv;

  if (pguess < SMALLP) { pguess = SMALLP; }

  /* Newton-Raphson iteration */
  int niter = 0;

  do {
    niter += 1;
    pold         = pguess;
    float fL     = fp(pguess, left, AL, BL, aL);
    float fR     = fp(pguess, right, AR, BR, aR);
    float dfpdpL = dfpdp(pguess, left, AL, BL, aL);
    float dfpdpR = dfpdp(pguess, right, AR, BR, aR);
    pguess       = pold - (fL + fR + delta_u) / (dfpdpL + dfpdpR);
    if (pguess < EPSILON_ITER) { pguess = SMALLP; }
    if (niter > 100) {
      printf(
        "Iteration for central pressure needs more than %d steps. "
        "Force-quitting iteration. Old-to-new ratio is %g\n",
        niter,
        fabsf(1.f - pguess / pold)
      );
      break;
    }
  } while (2.f * fabsf((pguess - pold) / (pguess + pold)) >= EPSILON_ITER);

  if (pguess <= SMALLP) { pguess = SMALLP; }

  *ustar = left->u[dim] - fp(pguess, left, AL, BL, aL);
  *pstar = pguess;

  /* debugmessage("p* found after %d iterations.", niter); */
  /* debugmessage("Got pstar = %12.8f, ustar = %12.8f", *pstar, *ustar); */
  /* debugmessage("rhoL %.3f uL %.3f pL %.3f rhoR %.3f uR %.3f pR %.3f dim
   * %d\n",  */
  /*     left->rho, left->u[dim], left->p, right->rho, right->u[dim], right->p,
   * dim); */
}

float fp(float pstar, pstate* s, float A, float B, float a) {
  /* -----------------------------------------------------*/
  /* Left/Right part of the pressure function             */
  /*------------------------------------------------------*/

  if (pstar > s->p) {
    /* we have a shock situation */
    return (pstar - s->p) * sqrtf(A / (pstar + B));
  }
  /* we have a rarefaction situation */
  return 2. * a / GM1 * (powf(pstar / s->p, BETA) - 1.);
}

float dfpdp(float pstar, pstate* s, float A, float B, float a) {
  /*---------------------------------------------------------------*/
  /* First derivative of Left/Right part of the pressure function  */
  /*---------------------------------------------------------------*/

  if (pstar > s->p) {
    /* we have a shock situation */
    return sqrtf(A / (pstar + B)) * (1. - 0.5 * (pstar - s->p) / (pstar + B));
  }
  /* we have a rarefaction situation */
  return 1. / (s->rho * a) * pow(pstar / s->p, -0.5 * GP1 / GAMMA);
}
