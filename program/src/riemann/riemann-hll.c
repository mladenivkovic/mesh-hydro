/* HLL Riemann Solver     */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#include "params.h"
#include "gas.h"
#include "riemann-hll.h"
#include "godunov.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern double gamma;
extern params pars;
extern double* x;
/* extern double t; */
extern cstate* flux;
extern cstate* u_old;
extern pstate* w_old;



/* =================================== */
void compute_fluxes(){
/* =================================== */
  /* compute the intercell fluxes      */
  /* fluxes are for conserved          */
  /* variables, not primitives!!!!     */
  /*-----------------------------------*/

  double SL = 0;
  double SR = 0;
  for (int i=NBC; i<pars.nx+NBCT; i++){
    cstate left = u_old[i-1];
    cstate right = u_old[i];
    compute_wave_speeds(left, right, &SL, &SR);
  
    pstate l = w_old[i-1];
    pstate r = w_old[i];
    cstate fl;
    cstate fr;
    fl.rho = l.rho * l.u;
    fl.rhou = fl.rho * l.u + l.p;
    fl.E = l.u*(energy(&l) + l.p);
    fr.rho = r.rho * r.u;
    fr.rhou = fr.rho * r.u + r.p;
    fr.E = r.u*(energy(&r) + r.p);

    if (SL >= 0){
      flux[i].rho = fl.rho;
      flux[i].rhou = fl.rhou;
      flux[i].E = fl.E;
    }
    else if (SL<=0 && SR>=0){
      flux[i].rho = Fhll(left.rho, right.rho, fl.rho, fr.rho, SL, SR);
      flux[i].rhou = Fhll(left.rhou, right.rhou, fl.rhou, fr.rhou, SL, SR);
      flux[i].E = Fhll(left.E, right.E, fl.E, fr.E, SL, SR);
    }
    else{
      flux[i].rho = fr.rho;
      flux[i].rhou = fr.rhou;
      flux[i].E = fr.E;
    }

    /* TODO: temporary checks */
    /* printf("Got flux for x=%lf: SL= %lf\t SR=%lf \t",x[i], SL, SR); */
    /* printf("%lf \t %lf \t %lf \n", flux[i].rho, flux[i].rhou , flux[i].E); */
    /* printf("%lf \t %lf \t %lf \n", u_old[pars.nx+NBC-1].rho, u_old[pars.nx+NBC-1].rhou , u_old[pars.nx+NBC-1].E); */
    /* printf("%lf \t %lf \t %lf \n", u_old[pars.nx+NBC].rho, u_old[pars.nx+NBC].rhou , u_old[pars.nx+NBC].E); */
    /* printf("%lf \t %lf \t %lf \n", w_old[pars.nx+NBC-1].rho, w_old[pars.nx+NBC-1].u , w_old[pars.nx+NBC-1].p); */
    /* printf("%lf \t %lf \t %lf \n", w_old[pars.nx+NBC].rho, w_old[pars.nx+NBC].u , w_old[pars.nx+NBC].p); */
    /* printf("-----------------------------\n"); */
  }


}




/* [> ============================================================================== <] */
/* void compute_uhll(cstate* left, cstate* right, cstate* u_hll, double SL, double SR){ */
/* [> ============================================================================== <] */
/*   [> Compute the HLL intermediate state                                           <] */
/*   [>------------------------------------------------------------------------------<] */
/*  */
/*   cstate FL, FR; */
/*   FL.rho = left->rho * left->u; */
/*   FL.rhou = left->rho * left->u * left->u + left->p; */
/*   FL.E = left->u * (energy(left) + left->p; */
/*   FR.rho = right->rho * right->u; */
/*   FR.rhou = right->rho * right->u * right->u + right->p; */
/*   FR.E = right->u * (energy(right) + right->p; */
/*    */
/*   u_hll->rho  = uhll(left->rho,  right->rho,  FL->rho,  FR->rho,  SL, SR); */
/*   u_hll->rhou = uhll(left->rhou, right->rhou, FL->rhou, FR->rhou, SL, SR); */
/*   u_hll->E    = uhll(left->E,    right->E,    FL->E,    FR->E,    SL, SR); */
/*  */
/* } */


/* ============================================================================= */
double Fhll(double ul, double ur, double fl, double fr, double SL, double SR){
/* ============================================================================= */
  /* compute hll Flux for a single component                                     */
  /*-----------------------------------------------------------------------------*/
  return (SR*fl - SL*fr + SL*SR*(ur - ul))/(SR - SL);
}


/* ========================================================================== */
void compute_wave_speeds(cstate left, cstate right, double* SL, double* SR){
/* ========================================================================== */
  /* Compute the Roe average wave speeds                                      */
  /*--------------------------------------------------------------------------*/

  double sqrtrl = sqrt(left.rho);
  double sqrtrr = sqrt(right.rho);

  double ul = left.rhou / left.rho;
  double ur = right.rhou / right.rho;

  double pl = (gamma-1)*(left.E - 0.5*left.rhou*ul);
  double pr = (gamma-1)*(right.E - 0.5*right.rhou*ur);

  double Hl = (left.E + pl)/left.rho;
  double Hr = (right.E + pr)/right.rho;

  double utilde = (sqrtrl * ul + sqrtrr*ur)/(sqrtrl+sqrtrr);
  double Htilde = (sqrtrl * Hl + sqrtrr*Hr)/(sqrtrl+sqrtrr);
  double temp = fabs(Htilde - 0.5*utilde*utilde); /* me trying to fix nans. Setting it to 0 also doesn't help.*/
  double atilde = sqrt((gamma-1)*(temp));

  *SL = utilde - atilde;
  *SR = utilde + atilde;

  /* TODO: temporary checks */
  /* printf("Computing wave speeds: utilde=%lf\tatilde=%lf\t rhol=%lf rhor=%lf\n", utilde, atilde, left.rho, right.rho); */

}
