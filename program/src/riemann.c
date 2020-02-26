/* functions that all Riemann solvers use */

/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */


#include <math.h>

#include "gas.h"
#include "params.h"
#include "utils.h"

extern params pars;



int riemann_has_vacuum(pstate *left, pstate *right, int dimension){
  /* ------------------------------------------------------------------------- 
   * Check whether we work with vacuum                    
   * returns true (1) if vacuum, 0 otherwise              
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * int dimension:   which fluid velocity dimension to use. 0: x, 1: y
   * ------------------------------------------------------------------------- */

  if (left->rho == 0){ return(1); }
  if (right->rho == 0){ return(1); }

  float delta_u = right->u[dimension] - left->u[dimension];
  float u_crit = 2./GM1 * (gas_soundspeed(left) + gas_soundspeed(right));

  if (delta_u < u_crit){
    return(0);
  } else {
    return(1);
  }
}






void riemann_compute_vacuum_solution(pstate* left, pstate* right, pstate* sol, 
      float xovert, float* wavevel, int dim){
  /* -------------------------------------------------------------------------
   * Solve the Riemann problem posed by a left and right state
   * and sample the solution at given x and t
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * pstate* sol:     pstate where solution will be written
   * float xovert:    x / t, point where solution shall be sampled
   * float* wavevel:  highest wave velocity from the problem
   * int dim:   which fluid velocity dim to use. 0: x, 1: y
   * ------------------------------------------------------------------------- */


  if (left->rho==0 && right->rho == 0){
    sol->rho = 0;
    sol->u[dim] = 0;
    sol->p = 0;
    return;
  }

  if (left->rho==0){
    /*------------------------*/
    /* Left vacuum state      */
    /*------------------------*/
    float aR = gas_soundspeed(right);
    float SR = right->u[dim] - 2*aR/GM1; /* vacuum front speed */
    *wavevel = fabs(SR);
    float SHR = right->u[dim] + aR;   /* speed of head of right rarefaction fan */
    if (fabs(SHR) > *wavevel) *wavevel = fabs(SHR);

    if (xovert <= SR){
      /* left vacuum */
      sol->rho = 0.;
      sol->u[dim] = SR;
      sol->p = 0.;
    } 
    else if (xovert < SHR){
      /* inside rarefaction */
      float precomp = pow(( 2. / GP1 - GM1OGP1 / aR *(right->u[dim] - xovert) ), (2/GM1));
      sol->rho = right->rho * precomp;
      sol->u[dim] = 2./ GP1 * (GM1HALF * right->u[dim] - aR + xovert);
      sol->p = right->p * pow(precomp, GAMMA);
    }
    else{
      /* original right pstate */
      sol->rho = right->rho;
      sol->u[dim] = right->u[dim];
      sol->p = right->p;
    }
  }

  else if (right->rho==0){
    /*------------------------*/
    /* Right vacuum state     */
    /*------------------------*/

    float aL = gas_soundspeed(left);
    float SL = left->u[dim] + 2*aL/GM1; /* vacuum front speed */
    *wavevel = fabs(SL);
    float SHL = left->u[dim] - aL;    /* speed of head of left rarefaction fan */
    if (fabs(SHL) > *wavevel) *wavevel = fabs(SHL);

    if (xovert >= SL){
      /* right vacuum */
      sol->rho = 0.;
      sol->u[dim] = SL;
      sol->p = 0.;
    }
    else if (xovert > SHL){
      /* inside rarefaction */
      float precomp = pow(( 2. / GP1 + GM1OGP1 / aL *(left->u[dim] - xovert) ), (2./GM1));
      sol->rho = left->rho * precomp;
      sol->u[dim] = 2./GP1 * (GM1HALF * left->u[dim] + aL + xovert);
      sol->p = left->p * pow(precomp, GAMMA);
    }
    else{
      /* original left pstate */
      sol->rho = left->rho;
      sol->u[dim] = left->u[dim];
      sol->p = left->p;
    }
  }
  else {
    /*------------------------*/
    /* Vacuum generating case */
    /*------------------------*/

    float aL = gas_soundspeed(left);
    float aR = gas_soundspeed(right);
    float SL = left->u[dim] + 2*aL/GM1; /* vacuum front speed */
    float SR = right->u[dim] - 2*aR/GM1; /* vacuum front speed */
    float SHL = left->u[dim] - aL;    /* speed of head of left rarefaction fan */
    float SHR = right->u[dim] + aR;   /* speed of head of right rarefaction fan */

    /* store max velocity */
    *wavevel = fabs(SL);
    if (fabs(SR) > *wavevel) *wavevel = fabs(SR);
    if (fabs(SHR) > *wavevel) *wavevel = fabs(SHR);
    if (fabs(SHL) > *wavevel) *wavevel = fabs(SHL);

    if (xovert <= SHL){
      /* left original pstate*/
      sol->rho = left->rho;
      sol->u[dim] = left->u[dim];
      sol->p = left->p;
    }
    else if (xovert < SL){
      /* inside rarefaction fan from right to left */
      float precomp = pow(( 2. / GP1 + GM1OGP1 / aL *(left->u[dim] - xovert) ), (2./GM1));
      sol->rho = left->rho * precomp;
      sol->u[dim] = 2./GP1 * (GM1HALF * left->u[dim] + aL + xovert);
      sol->p = left->p * pow(precomp, GAMMA);
    }
    else if (xovert < SR) {
      /* vacuum region */
      sol->rho = 0;
      sol->u[dim] = 0.5 * (SL + SR); /* just made something up here */
      sol->p = 0;
    }
    else if (xovert < SHR){
      /* inside rarefaction fan from left to right */
      float precomp = pow(( 2. / GP1 - GM1OGP1 / aR *(right->u[dim] - xovert) ), (2/GM1));
      sol->rho = right->rho * precomp;
      sol->u[dim] = 2./ GP1 * (GM1HALF * right->u[dim] - aR + xovert);
      sol->p = right->p * pow(precomp, GAMMA);
    }
    else{
      /* right original pstate */
      sol->rho = right->rho;
      sol->u[dim] = right->u[dim];
      sol->p = right->p;
    }
  }

  return;
}
