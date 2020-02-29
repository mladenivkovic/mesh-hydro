/* Two Rarefaction Approximate Riemann Solver */

/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */

#include "params.h"
#include "gas.h"
#include "riemann.h"
#include "utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern params pars;



void riemann_solve(pstate* left, pstate* right, pstate* sol, float xovert, int dimension){
  /* -------------------------------------------------------------------------
   * Solve the Riemann problem posed by a left and right state
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * pstate* sol:     pstate where solution will be written
   * float xovert:    x / t, point where solution shall be sampled
   * int dimension:   which fluid velocity dimension to use. 0: x, 1: y
   * ------------------------------------------------------------------------- */

  if (riemann_has_vacuum(left, right, dimension)){
    riemann_compute_vacuum_solution(left, right, sol, xovert, dimension);
  } else {
    float pstar = 0;
    float ustar = 0;
    riemann_compute_star_states(left, right, &pstar, &ustar, dimension);
    riemann_sample_solution(left, right, pstar, ustar, sol, xovert, dimension);
  }
}






void riemann_compute_star_states(pstate *left, pstate *right, float *pstar, float *ustar, int dim){
  /* ------------------------------------------------------------------------------------------
   * computes the star region pressure and velocity given the left and right pstates.         
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * float* pstar:    where pressure of star region will be written
   * float* ustar:    where velocity of star region will be written
   * int dimension:   which fluid velocity dimensionto use. 0: x, 1: y
   *------------------------------------------------------------------------------------------*/

  float aL = gas_soundspeed(left);
  float aLinv = 1./aL;
  float aR = gas_soundspeed(right);
  float aRinv = 1./aR;

  float pLRbeta = pow(left->p/right->p, BETA);

  *ustar = ( (pLRbeta  - 1.) / GM1HALF + left->u[dim] * aLinv * pLRbeta + right->u[dim]*aRinv ) /
              ( aRinv + aLinv * pLRbeta  );

  *pstar = 0.5 * (right->p * pow((1. + aRinv * GM1HALF * (*ustar - right->u[dim])), 1./BETA) +
                  left->p  * pow((1. + aLinv * GM1HALF * (left->u[dim]  - *ustar)), 1./BETA));

  debugmessage("Got pstar = %12.8f, ustar = %12.8f", *pstar, *ustar);
}







void riemann_sample_solution(pstate* left, pstate* right, float pstar, 
  float ustar, pstate* sol, float xovert, int dim){
  /*--------------------------------------------------------------------------------------------------
   * Compute the solution of the riemann problem at given time t and x, specified as xovert = x/t     
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * float pstar:     pressure of star region
   * float ustar:     velocity of star region
   * pstate* sol:     pstate where solution will be written
   * float xovert:    x / t, point where solution shall be sampled
   * int dim:         which fluid velocity direction to use. 0: x, 1: y
   *--------------------------------------------------------------------------------------------------*/

  if (xovert <= ustar){
    /*------------------------*/
    /* We're on the left side */
    /*------------------------*/
    float aL =  gas_soundspeed(left);
    float pstaroverpL = pstar/left->p;

    if (pstar <= left->p){
      /*------------------*/
      /* left rarefaction */
      /*------------------*/
      float SHL = left->u[dim] - aL;    /* speed of head of left rarefaction fan */
      if (xovert < SHL) {
        /* we're outside the rarefaction fan */
        sol->rho = left->rho;
        sol->u[dim] = left->u[dim];
        sol->p = left->p;
      }
      else {
        float astarL = aL * pow(pstaroverpL, BETA);
        float STL = ustar - astarL;  /* speed of tail of left rarefaction fan */
        if (xovert < STL){
          /* we're inside the fan */
          float precomp = pow(( 2. / GP1 + GM1OGP1 / aL *(left->u[dim] - xovert) ), (2./GM1));
          sol->rho = left->rho * precomp;
          sol->u[dim] = 2./GP1 * (GM1HALF * left->u[dim] + aL + xovert);
          sol->p = left->p * pow(precomp, GAMMA);
        }
        else{
          /* we're in the star region */
          sol->rho = left->rho*pow(pstaroverpL, ONEOVERGAMMA);
          sol->u[dim] = ustar;
          sol->p = pstar;
        }
      }
    }
    else{
      /*------------------*/
      /* left shock       */
      /*------------------*/
      float SL  = left->u[dim]  - aL * sqrtf(0.5 * GP1/GAMMA * pstaroverpL + BETA); /* left shock speed */
      if (xovert < SL){
        /* we're outside the shock */
        sol->rho = left->rho;
        sol->u[dim] = left->u[dim];
        sol->p = left->p;
      }
      else{
        /* we're in the star region */
        sol->rho = (pstaroverpL + GM1OGP1) / (GM1OGP1 * pstaroverpL + 1) * left->rho;
        sol->u[dim] = ustar;
        sol->p = pstar;
      }
    }
  }
  else{
    /*-------------------------*/
    /* We're on the right side */
    /*-------------------------*/
    float aR =  gas_soundspeed(right);
    float pstaroverpR = pstar / right->p;
    if (pstar <= right->p){

      /*-------------------*/
      /* right rarefaction */
      /*-------------------*/
      float SHR = right->u[dim] + aR;   /* speed of head of right rarefaction fan */
      if (xovert > SHR) {
        /* we're outside the rarefaction fan */
        sol->rho = right->rho;
        sol->u[dim] = right->u[dim];
        sol->p = right->p;
      }
      else {
        float astarR = aR * pow(pstaroverpR, BETA);
        float STR = ustar + astarR;  /* speed of tail of right rarefaction fan */
        if (xovert > STR){
          /* we're inside the fan */
          float precomp = pow(( 2. / GP1 - GM1OGP1 / aR *(right->u[dim] - xovert) ), (2/GM1));
          sol->rho = right->rho * precomp;
          sol->u[dim] = 2./ GP1 * (GM1HALF * right->u[dim] - aR + xovert);
          sol->p = right->p * pow(precomp, GAMMA);
        }
        else{
          /* we're in the star region */
          sol->rho = right->rho * pow(pstaroverpR, ONEOVERGAMMA);
          sol->u[dim] = ustar;
          sol->p = pstar;
        }
      }
    }
    else{
      /*------------------*/
      /* right shock      */
      /*------------------*/
      float SR  = right->u[dim] + aR*sqrtf(0.5*GP1/GAMMA * pstaroverpR + BETA); /* right shock speed */
      if (xovert > SR){
        /* we're outside the shock */
        sol->rho = right->rho;
        sol->u[dim] = right->u[dim];
        sol->p = right->p;
      }
      else{
        /* we're in the star region */
        sol->rho = (pstaroverpR + GM1OGP1) / (GM1OGP1 * pstaroverpR + 1) * right->rho;
        sol->u[dim] = ustar;
        sol->p = pstar;
      }
    }
  }

  return;
}
