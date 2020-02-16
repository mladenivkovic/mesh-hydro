/* Exact Riemann Solver */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#include "params.h"
#include "gas.h"
#include "riemann-exact.h"
#include "godunov.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern double gamma;
extern params pars;

/* extern double* x; */
extern pstate* w_old;
extern pstate* w_intercell;
extern cstate* flux;
extern cstate* u_old;
/* extern double t; */


/* =================================== */
void compute_fluxes(){
/* =================================== */
  /* compute the intercell fluxes      */
  /* fluxes are for conserved          */
  /* variables, not primitives!!!!     */
  /*-----------------------------------*/

  compute_intercell_states();

  for (int i=0; i<pars.nx+NBCT; i++){
    pstate ic = w_intercell[i];

    flux[i].rho = ic.rho * ic.u;
    flux[i].rhou = ic.rho * ic.u*ic.u + ic.p;
    flux[i].E = ic.u*(energy(&ic) + ic.p);
  }
}


/* ===================================== */
void compute_intercell_states(){
/* ===================================== */
  /* Compute intercell states, based on  */
  /* whether we have vacuum or not       */
  /* first compute primitives, then      */
  /* convert to conserved variables      */
  /*-------------------------------------*/

  pstate left, right, starL, starR;

  for(int i = 1; i<pars.nx+NBCT; i++){
    left = w_old[i-1];
    right = w_old[i];

    int vacuum = check_vacuum(&left, &right);

    if (vacuum){
      /* printf("calling vacuum for x=%lf\n", x[i]); */
      compute_riemann_vacuum(&left, &right, &w_intercell[i]);
    }
    else {
      compute_star_pstate(&left, &right, &starL, &starR);
      compute_riemann(&left, &right, &starL, &starR, &w_intercell[i]);
    }
  }
}





/* ====================================================== */
int check_vacuum(pstate *left, pstate *right){
/* ====================================================== */
  /* Check whether we work with vacuum                    */
  /* returns true (1) if vacuum, 0 otherwise              */
  /*------------------------------------------------------*/

  if (left->rho == 0 && left->p == 0){
    return(1);
  }
  if (right->rho == 0 && right->p == 0){
    return(1);
  }

  double delta_u = right->u - left->u;
  double u_crit = 2./(gamma - 1) * (soundspeed(left) + soundspeed(right));

  if (delta_u < u_crit){
    return(0);
  }
  else {
    return(1);
  }
}





/* ========================================================================================== */
void compute_star_pstate(pstate *left, pstate *right, pstate* starL, pstate* starR){
/* ========================================================================================== */
  /* computes the star pstate given the left and right pstates.                               */
  /*------------------------------------------------------------------------------------------*/

  double AL = 2. / ((gamma + 1) * left->rho);
  double AR = 2. / ((gamma + 1) * right->rho);
  double BL = (gamma - 1)/(gamma + 1) * left->p;
  double BR = (gamma - 1)/(gamma + 1) * right->p;
  double aL = soundspeed(left);
  double aR = soundspeed(right);

  double delta_u = right->u - left->u;

  double tolerance = 1e-6;

  double pguess, pold;
  double fL, fR, dfpdpL, dfpdpR;


  pguess = 0.5*(left->p + right->p);

  /* int niter = 0; */

  do {
    /* niter += 1; */
    fL = fp(pguess, left,  gamma, AL, BL, aL);
    fR = fp(pguess, right, gamma, AR, BR, aR);
    dfpdpL = dfpdp(pguess, left,  gamma, AL, BL, aL);
    dfpdpR = dfpdp(pguess, right, gamma, AR, BR, aR);
    pold = pguess;
    pguess -= fpfull(fL, fR, delta_u)/dfpdpfull(dfpdpL, dfpdpR);
    if (pguess<tolerance) pguess=tolerance;
    /* //printf("pstar iter %d pguess=%10.6lf, pold=%10.6lf\n", niter, pguess, pold); */
  }
  while (2*fabs((pguess-pold)/(pguess+pold)) >= tolerance);
  /* //printf("p* found after %d iterations.\n", niter); */


  starL->p = pguess;
  starR->p = pguess;

  fL = fp(pguess, left,  gamma, AL, BL, aL);
  fR = fp(pguess, right, gamma, AR, BR, aR);
  double u = 0.5*((left->u + right->u) + fR - fL);
  starL->u = u;
  starR->u = u;

  starL->rho = rho_star(left,  starL);
  starR->rho = rho_star(right, starR);
  
}




/* ==================================================== */
double rho_star(pstate *s, pstate *star){
/* ==================================================== */
  /* Compute the density in the star region             */
  /*----------------------------------------------------*/

  double pdiv = (star->p/s->p);
  if (star->p > s->p) {
    /* shocking matters */
    double gamfact = (gamma - 1)/(gamma + 1);
    return s->rho * (gamfact + pdiv) / (gamfact*pdiv + 1);
  } 
  else{
    /* rare occasions */
    return s->rho * pow(pdiv, (1./gamma));
  }
}





/* ========================================================== */
double fpfull(double fpL, double fpR, double delta_u){
/* ========================================================== */
  /* Full pressure function                                   */
  /*----------------------------------------------------------*/
  return fpL + fpR + delta_u;
}


/* ================================================== */
double dfpdpfull(double dfpdpL, double dfpdpR){
/* ================================================== */
  /* Full derivative of pressure function             */
  /*--------------------------------------------------*/
  return dfpdpL + dfpdpR;
}






/* =============================================================================== */
double fp(double pguess, pstate *s, double gamma, double A, double B, double a){
/* =============================================================================== */
  /* Left/Right part of the pressure function                                      */
  /*-------------------------------------------------------------------------------*/

  if (pguess > s->p){
    /* we have a shock situation */
    return (pguess - s->p)*sqrt(A/(pguess + B));
  }
  else{
    /* we have a rarefaction situation */
    return 2 * a / (gamma - 1) * ( pow(pguess/s->p, 0.5*(gamma-1)/gamma) - 1 );
  }
}


/* =============================================================================== */
double dfpdp(double pguess, pstate *s, double gamma, double A, double B, double a){
/* =============================================================================== */
  /* First derivative of Left/Right part of the pressure function                  */
  /*-------------------------------------------------------------------------------*/

  if (pguess > s->p){
    /* we have a shock situation */
    return sqrt(A/(pguess + B)) * (1 + 0.5 * (pguess - s->p)/(pguess + B));
  }
  else{
    /* we have a rarefaction situation */
    return 1./(s->rho * a) * ( pow(pguess/s->p, -0.5*(gamma+1)/gamma) );
  }
}








/* ================================================================================================== */
void compute_riemann(pstate* left, pstate* right, pstate* starL, pstate* starR, pstate* intercell){
/* ================================================================================================== */
  /* Compute the solution of the riemann problem at given time t for all x                            */
  /*--------------------------------------------------------------------------------------------------*/

  double S = 0;


  if (S <= starL->u){
    /*------------------------*/
    /* We're on the left side */
    /*------------------------*/
    double al =  soundspeed(left);
    if (starL->p <= left->p){
      /*------------------*/
      /* left rarefaction */
      /*------------------*/
      double SHL = left->u - al;    /* speed of head of left rarefaction fan */
      if (S < SHL) {
        /* we're outside the rarefaction fan */
        intercell->rho = left->rho;
        intercell->u = left->u;
        intercell->p = left->p;
      }
      else {
        double asl = soundspeed(starL);
        double STL = starL->u - asl;  /* speed of tail of left rarefaction fan */
        if (S < STL){
          /* we're inside the fan */
          intercell->rho = rho_fanL(left);
          intercell->u = u_fanL(left);
          intercell->p = p_fanL(left);
        }
        else{
          /* we're in the star region */
          intercell->rho = starL->rho;
          intercell->u = starL->u;
          intercell->p = starL->p;
        }
      }
    }
    else{
      /*------------------*/
      /* left shock       */
      /*------------------*/
      double SL  = left->u  - al*sqrt(0.5*((gamma+1) * starL->p/left->p + (gamma-1))/gamma); /* left shock speed */
      if (S<SL){
        /* we're outside the shock */
        intercell->rho = left->rho;
        intercell->u = left->u;
        intercell->p = left->p;
      }
      else{
        /* we're in the star region */
        intercell->rho = starL->rho;
        intercell->u = starL->u;
        intercell->p = starL->p;
      }
    }
  }
  else{
    /*-------------------------*/
    /* We're on the right side */
    /*-------------------------*/
    double ar =  soundspeed(right);
    if (starR->p <= right->p){
      /*-------------------*/
      /* right rarefaction */
      /*-------------------*/
      double SHR = right->u + ar;   /* speed of head of right rarefaction fan */
      if (S > SHR) {
        /* we're outside the rarefaction fan */
        intercell->rho = right->rho;
        intercell->u = right->u;
        intercell->p = right->p;
      }
      else {
        double asr = soundspeed(starR);
        double STR = starR->u + asr;  /* speed of tail of right rarefaction fan */
        if (S > STR){
          /* we're inside the fan */
          intercell->rho = rho_fanR(right);
          intercell->u = u_fanR(right);
          intercell->p = p_fanR(right);
        }
        else{
          /* we're in the star region */
          intercell->rho = starR->rho;
          intercell->u = starR->u;
          intercell->p = starR->p;
        }
      }
    }
    else{
      /*------------------*/
      /* right shock      */
      /*------------------*/
      double SR  = right->u + ar*sqrt(0.5*((gamma+1) * starR->p/right->p + (gamma-1))/gamma); /* right shock speed */
      if (S>SR){
        /* we're outside the shock */
        intercell->rho = right->rho;
        intercell->u = right->u;
        intercell->p = right->p;
      }
      else{
        /* we're in the star region */
        intercell->rho = starR->rho;
        intercell->u = starR->u;
        intercell->p = starR->p;
      }
    }
  }

  return;
}




/* ================================================================================ */
void compute_riemann_vacuum(pstate* left, pstate* right, pstate* intercell){
/* ================================================================================ */
  /* Compute the solution of the riemann problem at given time t for x = 0          */
  /*--------------------------------------------------------------------------------*/

  if (left->rho==0){
    /*------------------------*/
    /* Left vacuum state      */
    /*------------------------*/
    double ar = soundspeed(right);
    double SR = right->u - 2*ar/(gamma-1);
    double S = 0;

    if (S <= SR){
      /* left vacuum */
      intercell->rho = left->rho;
      intercell->u = SR;
      intercell->p = left->p;
    }
    else if (S < right->u + ar){
      /* inside rarefaction */
      intercell->rho = rho_fanR(right);
      intercell->u = u_fanR(right);
      intercell->p = p_fanR(right);
    }
    else{
      /* original right pstate */
      intercell->rho = right->rho;
      intercell->u = right->u;
      intercell->p = right->p;
    }

  }

  else if (right->rho==0){
    /*------------------------*/
    /* Right vacuum state     */
    /*------------------------*/

    double al = soundspeed(left);
    double SL = left->u + 2*al/(gamma-1);

    double S = 0;

    if (S >= SL){
      /* right vacuum */
      intercell->rho = right->rho;
      intercell->u = SL;
      intercell->p = right->p;
    }
    else if (S > left->u - al){
      /* inside rarefaction */
      intercell->rho = rho_fanL(left);
      intercell->u = u_fanL(left);
      intercell->p = p_fanL(left);
    }
    else{
      /* original left pstate */
      intercell->rho = left->rho;
      intercell->u = left->u;
      intercell->p = left->p;
    }

  }
  else {
    /*------------------------*/
    /* Vacuum generating case */
    /*------------------------*/

    double al = soundspeed(left);
    double ar = soundspeed(right);
    double SL = left->u + 2*al/(gamma-1);
    double SR = right->u - 2*ar/(gamma-1);

    for (int i=0; i<pars.nx; i++){

      double S = 0;

      if (S <= left->u-al){
        /* left original pstate*/
        intercell->rho = left->rho;
        intercell->u = left->u;
        intercell->p = left->p;
      }
      else if (S < SL){
        /* rarefaction fan from right to left */
        intercell->rho = rho_fanL(left);
        intercell->u = u_fanL(left);
        intercell->p = p_fanL(left);
      }
      else if (S < SR) {
        /* vacuum region */
        intercell->rho = 0;
        intercell->u = 0.5*(SL+SR); /* just made something up here */
        intercell->p = 0;
      }
      else if (S < right->u + ar){
        /* rarefaction fan from left to right */
        intercell->rho = rho_fanR(right);
        intercell->u = u_fanR(right);
        intercell->p = p_fanR(right);
      }
      else{
        /* right original pstate */
        intercell->rho = right->rho;
        intercell->u = right->u;
        intercell->p = right->p;
      }

    }
  }

  return;
}






/* ==================================== */
double rho_fanL(pstate* s){
/* ==================================== */
  /* Compute rho inside rarefaction fan */
  /* shortened version: S = x/t = 0     */
  /*------------------------------------*/
  return s->rho * pow((2. + (gamma-1)/soundspeed(s) * s->u)/(gamma+1), (2/(gamma-1)));
  /* return s->rho * pow((2./(gamma+1) + (gamma-1)/(gamma+1)/soundspeed(s) * (s->u - S)), (2/(gamma-1))); */
}

/* ==================================== */
double u_fanL(pstate* s){
/* ==================================== */
  /* Compute u inside rarefaction fan   */
  /* shortened version: S = x/t = 0     */
  /*------------------------------------*/
  return 2/(gamma+1) * (soundspeed(s) + 0.5*(gamma-1)*s->u);
  /* return 2/(gamma+1) * (soundspeed(s) + 0.5*(gamma-1)*s->u + S); */
}

/* ==================================== */
double p_fanL(pstate* s){
/* ==================================== */
  /* Compute p inside rarefaction fan   */
  /* shortened version: S = x/t = 0     */
  /*------------------------------------*/
  return s->p * pow((2 + (gamma-1)/soundspeed(s) * (s->u))/(gamma+1), (2*gamma/(gamma-1)));
  /* return s->p * pow((2/(gamma+1) + (gamma-1)/(gamma+1)/soundspeed(s) * (s->u - S)), (2*gamma/(gamma-1))); */
}





/* ==================================== */
double rho_fanR(pstate* s){
/* ==================================== */
  /* Compute rho inside rarefaction fan */
  /* shortened version: S = x/t = 0     */
  /*------------------------------------*/
  return s->rho * pow((2./(gamma+1) - (gamma-1)/(gamma+1)/soundspeed(s) * s->u), (2./(gamma-1)));
  /* return s->rho * pow((2./(gamma+1) - (gamma-1)/(gamma+1)/soundspeed(s) * (s->u - S)), (2./(gamma-1))); */
}

/* ==================================== */
double u_fanR(pstate* s){
/* ==================================== */
  /* Compute u inside rarefaction fan   */
  /* shortened version: S = x/t = 0     */
  /*------------------------------------*/
  return 2/(gamma+1) * (-soundspeed(s) + 0.5*(gamma-1)*s->u);
  /* return 2/(gamma+1) * (-soundspeed(s) + 0.5*(gamma-1)*s->u + S); */
}

/* ==================================== */
double p_fanR(pstate* s){
/* ==================================== */
  /* Compute p inside rarefaction fan   */
  /* shortened version: S = x/t = 0     */
  /*------------------------------------*/
  return s->p * pow((2./(gamma+1) - (gamma-1)/(gamma+1)/soundspeed(s) * (s->u)), (2*gamma/(gamma-1)));
  /* return s->p * pow((2./(gamma+1) - (gamma-1)/(gamma+1)/soundspeed(s) * (s->u - S)), (2*gamma/(gamma-1))); */
}
