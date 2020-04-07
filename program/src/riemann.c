/* functions that all Riemann solvers use */

/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */


#include <math.h>
#include <stddef.h>

#include "defines.h"
#include "gas.h"
#include "params.h"
#include "riemann.h"
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

  if (left->rho <= SMALLRHO){ return(1); }
  if (right->rho <= SMALLRHO){ return(1); }

  float delta_u = right->u[dimension] - left->u[dimension];
  float u_crit = 2./GM1 * (gas_soundspeed(left) + gas_soundspeed(right));

  if (delta_u < u_crit){
    return(0);
  } else {
    return(1);
  }
}






void riemann_compute_vacuum_solution(pstate* left, pstate* right, pstate* sol, 
      float xovert, int dim){
  /* -------------------------------------------------------------------------
   * Solve the Riemann problem posed by a left and right state
   * and sample the solution at given x and t
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * pstate* sol:     pstate where solution will be written
   * float xovert:    x / t, point where solution shall be sampled
   * int dim:   which fluid velocity dim to use. 0: x, 1: y
   * ------------------------------------------------------------------------- */


  if (left->rho <= SMALLRHO && right->rho == SMALLRHO){
    sol->rho = SMALLRHO;
    sol->u[dim] = SMALLU;
    sol->p = SMALLP;
    return;
  }

  if (left->rho <= SMALLRHO){
    /*------------------------*/
    /* Left vacuum state      */
    /*------------------------*/
    float aR = gas_soundspeed(right);
    float SR = right->u[dim] - 2*aR/GM1; /* vacuum front speed */
    float SHR = right->u[dim] + aR;   /* speed of head of right rarefaction fan */

    if (xovert <= SR){
      /* left vacuum */
      sol->rho = SMALLRHO;
#ifdef USE_AS_RIEMANN_SOLVER
      sol->u[dim] = SR;
#else
      sol->u[dim] = SMALLU;
#endif
      sol->p = SMALLP;
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

  else if (right->rho<=SMALLRHO){
    /*------------------------*/
    /* Right vacuum state     */
    /*------------------------*/

    float aL = gas_soundspeed(left);
    float SL = left->u[dim] + 2*aL/GM1; /* vacuum front speed */
    float SHL = left->u[dim] - aL;    /* speed of head of left rarefaction fan */

    if (xovert >= SL){
      /* right vacuum */
      sol->rho = SMALLRHO;
#ifdef USE_AS_RIEMANN_SOLVER
      sol->u[dim] = SL;
#else
      sol->u[dim] = SMALLU;
#endif
      sol->p = SMALLP;
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
      sol->rho = SMALLRHO;
#ifdef USE_AS_RIEMANN_SOLVER
      sol->u[dim] = 0.5*(SL + SR);
#else
      sol->u[dim] = SMALLU;
#endif
      sol->p = SMALLP;
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






void riemann_get_full_solution(pstate* left, pstate* right, float* S, cstate* fluxes, float* delta_q, int dim){
  /*-----------------------------------------------------------------------------------------------------------------
   *
   * TODO: DOX
   * ---------------------------------------------------------------------------------------------------------------- */

  /* First, get the star velocity and pressure */
  float pstar, ustar;
  riemann_compute_star_states(left, right, &pstar, &ustar, dim);

  /* Now fill them into left and right star state variables */
  pstate star_left;
  gas_init_pstate(&star_left);
  star_left.u[dim] = ustar;
  star_left.u[(dim+1) % 2] = left->u[(dim+1) % 2];
  star_left.p = pstar;
  float pstaroverpL = pstar / left->p;
  if (pstar <= left->p){
    /* rarefaction star state */
    star_left.rho = left->rho * pow(pstaroverpL, ONEOVERGAMMA);
  } else {
    /* shock star state */
    star_left.rho = (pstaroverpL + GM1OGP1) / (GM1OGP1 * pstaroverpL + 1.) * left->rho;
  }

  pstate star_right;
  gas_init_pstate(&star_right);
  star_right.u[dim] = ustar;
  star_right.u[(dim+1) % 2] = right->u[(dim+1) % 2];
  star_right.p = pstar;
  float pstaroverpR = pstar / right->p;
  if (pstar <= right->p){
    /* rarefaction star state */
    star_right.rho = right->rho * pow(pstaroverpR, ONEOVERGAMMA);
  } else{
    /* shock star state */
    star_right.rho = (pstaroverpR + GM1OGP1) / (GM1OGP1 * pstaroverpR + 1.) * right->rho;
  }





  /* get wave speeds. We always expect 3 waves. */

  float aL =  gas_soundspeed(left);
  if (pstar <= left->p){
    /* left rarefaction head velocity */
    S[0] = (left->u[dim] - aL);
  } else {
    /* left shock velocity */
    S[0] = (left->u[dim] - aL * sqrtf(0.5 * GP1/GAMMA * pstaroverpL + BETA));
  }

  S[1] = ustar;

  float aR =  gas_soundspeed(right);
  if (pstar <= right->p){
    /* right rarefaction head velocity */
    S[2] = (right->u[dim] + aR);
  } else {
    /* right shock velocity */
    S[2] = (right->u[dim] + aR * sqrtf(0.5 * GP1/GAMMA * pstaroverpR + BETA));
  }



  /* find solution at x = 0 in case we need it */
  pstate state_at_x_zero;
  gas_init_pstate(&state_at_x_zero);
  riemann_sample_solution(left, right, pstar, ustar, &state_at_x_zero, /*xovert =*/0., dim);

  /* [> find which of the four states is at x = 0. We need it in case we have sonic rarefactions. <] */
  /* pstate* state_at_x_zero = NULL; */
  /*  */
  /* [> printf("Velocities are %f %f %f ; state at x=0 is ", Sk[0], Sk[1], Sk[2]); <] */
  /* if (S[0] > 0.){ */
  /*   state_at_x_zero = left; */
  /*   [> printf("left\n"); <] */
  /* } else { */
  /*   [> left wave is negative, so find first positive velocity <] */
  /*   if (S[1] >= 0){ */
  /*     state_at_x_zero = &star_left; */
  /*     [> printf("star left\n"); <] */
  /*   } else if (S[2] >= 0){ */
  /*     state_at_x_zero = &star_right; */
  /*     [> printf("star right\n"); <] */
  /*   } else { */
  /*     state_at_x_zero = right; */
  /*     [> printf("right\n"); <] */
  /*   } */
  /* } */



  /* get fluxes from the states. For 1D Euler equations, there will always be 4 fluxes to consider. */

  /* state 1 */
  gas_get_cflux_from_pstate(left,  &fluxes[0], dim);

  /* TODO: document that we're putting rarefactions together with neighbour states */

  /* state 2 */
  if (pstar <= left->p){
    /* rarefaction. But is it sonic? */
    float astarL = gas_soundspeed(&star_left);
    float tailL = ustar - astarL;
    if (tailL*S[0] < 0.){
      /* if head and tail speeds have different sign, we have asonic rarefaction. */
      gas_get_cflux_from_pstate(&state_at_x_zero, &fluxes[1], dim);
    } else {
      /* Non-sonic rarefaction */
      gas_get_cflux_from_pstate(&star_left, &fluxes[1], dim);
    }
  } else {
    /* Shock */
    gas_get_cflux_from_pstate(&star_left, &fluxes[1], dim);
  }
  
  /* state 3 */
  if (pstar <= right->p){
    /* rarefaction. But is it sonic? */
    float astarR = gas_soundspeed(&star_right);
    float tailR = ustar + astarR;
    if (tailR*S[2] < 0.){
      /* if head and tail speeds have different sign, we have asonic rarefaction. */
      gas_get_cflux_from_pstate(&state_at_x_zero, &fluxes[2], dim);
    } else {
      /* Non-sonic rarefaction */
      gas_get_cflux_from_pstate(&star_right, &fluxes[2], dim);
    }
  } else {
    /* shock */
    gas_get_cflux_from_pstate(&star_right, &fluxes[2], dim);
  }

  /* state 4 */
  gas_get_cflux_from_pstate(right, &fluxes[3], dim);


  delta_q[0] = star_left.rho - left->rho;
  delta_q[1] = star_right.rho - star_left.rho;
  delta_q[2] = right->rho - star_right.rho;
}
