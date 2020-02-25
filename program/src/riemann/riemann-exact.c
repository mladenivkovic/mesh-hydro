/* Exact Riemann Solver */

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



/* ============================================================================== */
void riemann_solve(pstate* left, pstate* right, pstate* sol, 
      float xovert, float* wavevel, int dimension){
/* ============================================================================== */
  /* Solve the Riemann problem posed by a left and right state
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * pstate* sol:     pstate where solution will be written
   * float xovert:    x / t, point where solution shall be sampled
   * float* wavevel:  highest wave velocity from the problem
   * int dimension:   which fluid velocity dimension to use. 0: x, 1: y
   * ------------------------------------------------------------------------- */

    int vacuum = check_vacuum(&left, &right);

    if (vacuum){
      compute_riemann_vacuum(&left, &right, &w_intercell[i]);
    }
    else {
      float pstar = 0;
      float ustar = 0;
      riemann_compute_star_states(left, right, &pstar, &ustar, dimension);
      riemann_sample_solution(left, right, pstar, ustar, sol, xovert, wavevel, dimension);
    }
}






/* ========================================================================================== */
void riemann_compute_star_states(pstate *left, pstate *right, 
      float *pstar, float *ustar, int dim){
/* ========================================================================================== */
  /* computes the star region pressure and velocity given the left and right pstates.         
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
  float aL = gas_soundspeed(left->p, left->rho);
  float aR = gas_soundspeed(right->p, right->rho);

  float delta_u = right->u[dim] - left->u[dim];

  float pguess, pold;


  /* Find initial guess for star pressure */
  float ppv = 0.5*(left->p + right->p) - 0.125 * delta_u *
        (left->rho + right->rho) * (aL + aR);
  pguess = ppv;

  if (pguess < EPSILON_ITER) pguess = EPSILON_ITER;

  /* Newton-Raphson iteration */
  int niter = 0;

  do {
    niter += 1;
    pold = pguess;
    float fL = fp(pguess, left,  AL, BL, aL);
    float fR = fp(pguess, right, AR, BR, aR);
    float dfpdpL = dfpdp(pguess, left,  AL, BL, aL);
    float dfpdpR = dfpdp(pguess, right, AR, BR, aR);
    pguess = pold - (fL + fR + delta_u)/(dfpdpL + dfpdpR);
    if (pguess < EPSILON_ITER) pguess = EPSILON_ITER;
  }
  while (2*fabs((pguess-pold)/(pguess+pold)) >= EPSILON_ITER);

  *ustar = left->u[dim] - fp(pguess, left,  AL, BL, aL);
  *pstar = pguess;

  debugmessage("p* found after %d iterations.\n", niter);
  debugmessage("Got pstar = %12.8f, ustar = %12.8f\n", *pstar, *ustar);
  
}









/* =============================================================================== */
float fp(float pstar, pstate *s, float A, float B, float a){
/* =============================================================================== */
  /* Left/Right part of the pressure function                                      */
  /*-------------------------------------------------------------------------------*/

  if (pstar > s->p){
    /* we have a shock situation */
    return (pstar - s->p)*sqrtf(A/(pstar + B));
  }
  else{
    /* we have a rarefaction situation */
    return 2 * a / GM1 * (pow(pstar/s->p, ALPHA) - 1);
  }
}





/* =============================================================================== */
float dfpdp(float pstar, pstate *s, float A, float B, float a){
/* =============================================================================== */
  /* First derivative of Left/Right part of the pressure function                  */
  /*-------------------------------------------------------------------------------*/

  if (pstar > s->p){
    /* we have a shock situation */
    return sqrtf(A/(pstar + B)) * (1. - 0.5 * (pstar - s->p)/(pstar + B));
  }
  else{
    /* we have a rarefaction situation */
    return 1./(s->rho * a) * pow(pstar/s->p, -0.5*GP1/GAMMA);
  }
}








/* ================================================================================================== */
void riemann_sample_solution(pstate* left, pstate* right, 
    float pstar, float ustar, pstate* sol, float xovert, float* wavevel, int dim){
/* ================================================================================================== */
  /* Compute the solution of the riemann problem at given time t and x, specified as xovert = x/t     
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * float pstar:     pressure of star region
   * float ustar:     velocity of star region
   * pstate* sol:     pstate where solution will be written
   * float xovert:    x / t, point where solution shall be sampled
   * float* wavevel:  highest wave velocity from the problem
   * int dim:         which fluid velocity direction to use. 0: x, 1: y
   *--------------------------------------------------------------------------------------------------*/

  if (xovert <= ustar){
    /*------------------------*/
    /* We're on the left side */
    /*------------------------*/
    float aL =  gas_soundspeed(left->p, left->rho);
    float pstaroverpL = pstar/left->p;

    if (pstar <= left->p){
      /*------------------*/
      /* left rarefaction */
      /*------------------*/
      float SHL = left->u[dim] - aL;    /* speed of head of left rarefaction fan */
      *wavevel = fabs(SHL);
      if (xovert < SHL) {
        /* we're outside the rarefaction fan */
        sol->rho = left->rho;
        sol->u[dim] = left->u[dim];
        sol->p = left->p;
      }
      else {
        float astarL = aL * pow(pstaroverpL, ALPHA);
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
      float SL  = left->u[dim]  - aL * sqrtf(0.5 * GP1/GAMMA * pstaroverpL + ALPHA); /* left shock speed */
      *wavevel = fabs(SL);
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
    float aR =  gas_soundspeed(right->p, right->rho);
    float pstaroverpR = pstar / right->p;
    if (pstar <= right->p){

      /*-------------------*/
      /* right rarefaction */
      /*-------------------*/
      float SHR = right->u[dim] + aR;   /* speed of head of right rarefaction fan */
      *wavevel = fabs(SHR);
      if (xovert > SHR) {
        /* we're outside the rarefaction fan */
        sol->rho = right->rho;
        sol->u[dim] = right->u[dim];
        sol->p = right->p;
      }
      else {
        float astarR = aR * pow(pstaroverpR, ALPHA);
        float STR = ustar + astarR;  /* speed of tail of right rarefaction fan */
        if (xovert > STR){
          /* we're inside the fan */
          float precomp = pow(( 2. / GP1 - GM1OGP1 / aR *(right->u[dim] - xovert) ), (2/GM1));
          sol->rho = left->rho * precomp;
          sol->u[dim] = 2./ GP1 * (GM1HALF * right->u[dim] - aR + xovert);
          sol->p = left->p * pow(precomp, GAMMA);
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
      float SR  = right->u[dim] + aR*sqrtf(0.5*GP1/GAMMA * pstaroverpR + ALPHA); /* right shock speed */
      *wavevel = fabs(SR);
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




/* [> ================================================================================ <] */
/* void compute_riemann_vacuum(pstate* left, pstate* right, pstate* intercell){ */
/* [> ================================================================================ <] */
/*   [> Compute the solution of the riemann problem at given time t for x = 0          <] */
/*   [>--------------------------------------------------------------------------------<] */
/*  */
/*   if (left->rho==0){ */
/*     [>------------------------<] */
/*     [> Left vacuum state      <] */
/*     [>------------------------<] */
/*     float ar = gas_soundspeed(right); */
/*     float SR = right->ux - 2*ar/(gamma-1); */
/*     float S = 0; */
/*  */
/*     if (S <= SR){ */
/*       [> left vacuum <] */
/*       intercell->rho = left->rho; */
/*       intercell->ux = SR; */
/*       intercell->p = left->p; */
/*     } */
/*     else if (S < right->ux + ar){ */
/*       [> inside rarefaction <] */
/*       intercell->rho = rho_fanR(right); */
/*       intercell->ux = u_fanR(right); */
/*       intercell->p = p_fanR(right); */
/*     } */
/*     else{ */
/*       [> original right pstate <] */
/*       intercell->rho = right->rho; */
/*       intercell->ux = right->ux; */
/*       intercell->p = right->p; */
/*     } */
/*  */
/*   } */
/*  */
/*   else if (right->rho==0){ */
/*     [>------------------------<] */
/*     [> Right vacuum state     <] */
/*     [>------------------------<] */
/*  */
/*     float al = gas_soundspeed(left); */
/*     float SL = left->ux + 2*al/(gamma-1); */
/*  */
/*     float S = 0; */
/*  */
/*     if (S >= SL){ */
/*       [> right vacuum <] */
/*       intercell->rho = right->rho; */
/*       intercell->ux = SL; */
/*       intercell->p = right->p; */
/*     } */
/*     else if (S > left->ux - al){ */
/*       [> inside rarefaction <] */
/*       intercell->rho = rho_fanL(left); */
/*       intercell->ux = u_fanL(left); */
/*       intercell->p = p_fanL(left); */
/*     } */
/*     else{ */
/*       [> original left pstate <] */
/*       intercell->rho = left->rho; */
/*       intercell->ux = left->ux; */
/*       intercell->p = left->p; */
/*     } */
/*  */
/*   } */
/*   else { */
/*     [>------------------------<] */
/*     [> Vacuum generating case <] */
/*     [>------------------------<] */
/*  */
/*     float al = gas_soundspeed(left); */
/*     float ar = gas_soundspeed(right); */
/*     float SL = left->ux + 2*al/(gamma-1); */
/*     float SR = right->ux - 2*ar/(gamma-1); */
/*  */
/*     for (int i=0; i<pars.nx; i++){ */
/*  */
/*       float S = 0; */
/*  */
/*       if (S <= left->ux-al){ */
/*         [> left original pstate<] */
/*         intercell->rho = left->rho; */
/*         intercell->ux = left->ux; */
/*         intercell->p = left->p; */
/*       } */
/*       else if (S < SL){ */
/*         [> rarefaction fan from right to left <] */
/*         intercell->rho = rho_fanL(left); */
/*         intercell->ux = u_fanL(left); */
/*         intercell->p = p_fanL(left); */
/*       } */
/*       else if (S < SR) { */
/*         [> vacuum region <] */
/*         intercell->rho = 0; */
/*         intercell->ux = 0.5*(SL+SR); [> just made something up here <] */
/*         intercell->p = 0; */
/*       } */
/*       else if (S < right->ux + ar){ */
/*         [> rarefaction fan from left to right <] */
/*         intercell->rho = rho_fanR(right); */
/*         intercell->ux = u_fanR(right); */
/*         intercell->p = p_fanR(right); */
/*       } */
/*       else{ */
/*         [> right original pstate <] */
/*         intercell->rho = right->rho; */
/*         intercell->ux = right->ux; */
/*         intercell->p = right->p; */
/*       } */
/*  */
/*     } */
/*   } */
/*  */
/*   return; */
/* } */

