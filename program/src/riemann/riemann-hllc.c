/* HLLC Riemann Solver */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#include "gas.h"
#include "params.h"
#include "riemann.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void riemann_solve_hllc(pstate* left, pstate* right, cstate* sol, float xovert, int dimension){
  /* -------------------------------------------------------------------------
   * Solve the Riemann problem posed by a left and right state
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * cstate* sol:     cstate where solution (conserved FLUX) will be written
   * float xovert:    x / t, point where solution shall be sampled
   * int dimension:   which fluid velocity dimension to use. 0: x, 1: y
   * ------------------------------------------------------------------------- */

  if (riemann_has_vacuum(left, right, dimension)){
    /* the vacuum solver wants a pstate as the argument for the solution.
     * so give him one, and later translate it back to the conserved flux. */
    pstate pstate_sol;
    gas_init_pstate(&pstate_sol);
    riemann_compute_vacuum_solution(left, right, &pstate_sol, xovert, dimension);
    gas_get_cflux_from_pstate(&pstate_sol, sol, dimension);

  } else {

    float SL = 0;
    float SR = 0;
    riemann_compute_wave_speed_estimates(left, right, &SL, &SR, dimension);
    riemann_sample_solution(left, right, SL, SR, sol, xovert, dimension);
  }
}






void riemann_compute_wave_speed_estimates(pstate* left, pstate* right, float* SL, float* SR, int dimension){
  /*--------------------------------------------------------------------------------------------------------
   * Get estimates for the left and right HLLC wave speed.
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * int dimension:   which fluid velocity dimension to use. 0: x, 1: y
   * ------------------------------------------------------------------------------------------------------- */

  float aL = gas_soundspeed(left);
  float aR = gas_soundspeed(right);

  float delta_u = right->u[dimension] - left->u[dimension];
  float pstar = 0.5*(left->p + right->p) - 0.125 * delta_u *
        (left->rho + right->rho) * (aL + aR);
  if (pstar < 0) pstar = SMALLP;


  *SL = left->u[dimension] - aL * qLR(pstar, left->p);
  *SR = right->u[dimension] + aR * qLR(pstar, right->p);

}





float qLR(float pstar, float pLR){
  /*--------------------------------------------------
   * Compute q_{L,R} needed for the wave speed
   * estimate.
   * pstar:   (estimated) pressure of the star state
   * pLR:     left or right pressure, depending whether
   *          you want q_L or q_R
   * ------------------------------------------------- */

  if (pstar > pLR){
    /* shock relation */
    return (sqrtf(1. + BETA * (pstar/pLR - 1)));
  } else{
    /* rarefaction relation */
    return(1.); 
  }
}






void riemann_sample_solution(pstate* left, pstate* right, float SL, 
  float SR, cstate* sol, float xovert, int dim){
  /*--------------------------------------------------------------------------------------------------
   * Compute the solution of the riemann problem at given time t and x, specified as xovert = x/t     
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * float pstar:     pressure of star region
   * float ustar:     velocity of star region
   * cstate* sol:     cstate where solution (conserved FLUX) will be written
   * float xovert:    x / t, point where solution shall be sampled
   * int dim:         which fluid velocity direction to use. 0: x, 1: y
   *--------------------------------------------------------------------------------------------------*/

  float SLMUL = SL - left->u[dim];
  float SRMUR = SR - right->u[dim];
  float Sstar = ( right->p - left->p + left->rho * left->u[dim] * SLMUL - right->rho * right->u[dim] * SRMUR) /
      ( left->rho * SLMUL - right->rho * SRMUR );


  /* --------------------------------------- */ 
  /* compute left and right conserved states */
  /* --------------------------------------- */ 
  cstate UL;
  gas_init_cstate(&UL);
  gas_prim_to_cons(left, &UL);

  cstate UR;
  gas_init_cstate(&UR);
  gas_prim_to_cons(right, &UR);





  /* -------------------------------------------- */ 
  /* compute left and right conserved star states */
  /* -------------------------------------------- */ 

  cstate UstarL;
  gas_init_cstate(&UstarL);

  float lcomp = left->rho * SLMUL / (SL - Sstar);
  UstarL.rho = lcomp;
  UstarL.rhou[dim] = lcomp * Sstar;
  UstarL.rhou[(dim + 1) % 2] = lcomp * left->u[(dim + 1) % 2];

  if (left->rho > 0) {
    float EL= 0.5 * left->rho * (left->u[0] * left->u[0] + left->u[1] * left->u[1]) + left->p / GM1;
    UstarL.E = lcomp * ( (EL / left->rho ) + (Sstar - left->u[dim]) * 
        (Sstar + left->p / (left->rho * SLMUL)));
  } else {
    UstarL.E = SMALLRHO * SMALLU * SMALLU;
  }


  cstate UstarR;
  gas_init_cstate(&UstarR);

  float rcomp = right->rho * SRMUR / (SR - Sstar);
  UstarR.rho = rcomp;
  UstarR.rhou[dim] = rcomp * Sstar;
  UstarR.rhou[(dim + 1) % 2] = rcomp * right->u[(dim + 1) % 2];

  if (right->rho > 0) {
    float ER= 0.5 * right->rho * (right->u[0] * right->u[0] + right->u[1] * right->u[1]) + right->p / GM1;
    UstarR.E = rcomp * ( (ER / right->rho ) + (Sstar - right->u[dim]) * 
        (Sstar + right->p / (right->rho * SRMUR)));
  } else {
    UstarR.E = SMALLRHO * SMALLU * SMALLU;
  }





  /* ----------------------------------- */
  /* Compute left and right fluxes       */
  /* ----------------------------------- */
  cstate FL;
  gas_init_cstate(&FL);
  gas_get_cflux_from_cstate(&UL, &FL, dim);

  cstate FR;
  gas_init_cstate(&FR);
  gas_get_cflux_from_cstate(&UR, &FR, dim);







  /* ----------------------------------- */
  /* Compute left and right star fluxes  */
  /* ----------------------------------- */

  cstate FstarL;
  gas_init_cstate(&FstarL);

  FstarL.rho = FL.rho + SL * (UstarL.rho - UL.rho);
  FstarL.rhou[0] = FL.rhou[0] + SL * (UstarL.rhou[0] - UL.rhou[0]);
  FstarL.rhou[1] = FL.rhou[1] + SL * (UstarL.rhou[1] - UL.rhou[1]);
  FstarL.E = FL.E + SL * (UstarL.E - UL.E);


  cstate FstarR;
  gas_init_cstate(&FstarR);

  FstarR.rho = FR.rho + SL * (UstarR.rho - UR.rho);
  FstarR.rhou[0] = FR.rhou[0] + SL * (UstarR.rhou[0] - UR.rhou[0]);
  FstarR.rhou[1] = FR.rhou[1] + SL * (UstarR.rhou[1] - UR.rhou[1]);
  FstarR.E = FR.E + SL * (UstarR.E - UR.E);





  /* ---------------------------- */
  /* finally, sample the solution */
  /* ---------------------------- */

  if (xovert <= SL){
    /* solution is F_L */
    sol->rho = FL.rho;
    sol->rhou[0] = FL.rhou[0];
    sol->rhou[1] = FL.rhou[1];
    sol->E = FL.E;
  } else if (xovert <= Sstar){
    /* solution is F*_L */
    sol->rho = FstarL.rho;
    sol->rhou[0] = FstarL.rhou[0];
    sol->rhou[1] = FstarL.rhou[1];
    sol->E = FstarL.E;
  } else if ( xovert <= SR ){
    /* solution is F*_R */
    sol->rho = FstarR.rho;
    sol->rhou[0] = FstarR.rhou[0];
    sol->rhou[1] = FstarR.rhou[1];
    sol->E = FstarR.E;
  } else {
    /* solution is F_R */
    sol->rho = FR.rho;
    sol->rhou[0] = FR.rhou[0];
    sol->rhou[1] = FR.rhou[1];
    sol->E = FR.E;
  }


  return;
}
