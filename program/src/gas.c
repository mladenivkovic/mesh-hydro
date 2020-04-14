/* ideal gas related stuff */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "gas.h"
#include "defines.h"



void gas_init_pstate(pstate* s){
  /*-------------------------------------------------*/
  /* This function sets the pstate to zero           */
  /*-------------------------------------------------*/

  s->rho  = 0;
  s->u[0] = 0;
  s->u[1] = 0;
  s->p    = 0;
} 




void gas_init_cstate(cstate *s){
  /*-------------------------------------------------*/
  /* This function sets the cstate to zero           */
  /*-------------------------------------------------*/

  s->rho   = 0;
  s->rhou[0] = 0;
  s->rhou[1] = 0;
  s->E     = 0;
}




void gas_prim_to_cons(pstate* p, cstate* c){
  /* --------------------------------------------------------
   * Compute the conserved state vector of a given
   * primitive state
   * -------------------------------------------------------- */

  c->rho = p->rho;
  c->rhou[0] = p->rho * p->u[0];
  c->rhou[1] = p->rho * p->u[1];
  c->E = 0.5 * p->rho * (p->u[0] * p->u[0] + p->u[1] * p->u[1]) + p->p / GM1;
}




void gas_cons_to_prim(cstate *c, pstate* p){
  /* --------------------------------------------------------
   * Compute the primitive state vector of a given
   * conserved state
   * -------------------------------------------------------- */

  if (c->rho <= SMALLRHO){
    /* exception handling for vacuum */
    p->rho = SMALLRHO;
    p->u[0] = SMALLU;
    p->u[1] = SMALLU;
    p->p = SMALLP;
  } else {
    p->rho = c->rho;
    p->u[0] = c->rhou[0]/c->rho;
    p->u[1] = c->rhou[1]/c->rho;
    p->p = GM1 * ( c->E - 0.5 * (c->rhou[0] * c->rhou[0] + c->rhou[1] * c->rhou[1] ) / c->rho );
    /* do some exception handling. Sometimes the time step is too large, and we end up with negative pressures. */
    if (p->p <= SMALLP) p->p = SMALLP;
  }
}



void gas_get_cflux_from_pstate(pstate *p, cstate *f, int dimension){
  /* -----------------------------------------------------------
   * Compute the flux of conserved variables of the Euler
   * equations given a primitive state vector 
   *
   * The flux is not an entire tensor for 3D Euler equations, but
   * correpsonds to the dimensionally split vectors F, G as 
   * described in the "Euler equations in 2D" section of the
   * documentation TeX files.
   * That's why you need to specify the dimension.
   * ----------------------------------------------------------- */

  f->rho = p->rho * p->u[dimension];
  f->rhou[dimension] = p->rho * p->u[dimension] * p->u[dimension] + p->p;
  f->rhou[(dimension+1) % 2] = p->rho * p->u[0] * p->u[1];
  float E = 0.5 * p->rho * (p->u[0] * p->u[0] + p->u[1] * p->u[1]) + p->p / GM1;
  f->E = (E + p->p) * p->u[dimension];
}




void gas_get_cflux_from_cstate(cstate *c, cstate *f, int dimension){
  /* -----------------------------------------------------------
   * Compute the flux of conserved variables of the Euler
   * equations given a conserved state vector 
   *
   * The flux is not an entire tensor for 3D Euler equations, but
   * correpsonds to the dimensionally split vectors F, G as 
   * described in the "Euler equations in 2D" section of the
   * documentation TeX files.
   * That's why you need to specify the dimension.
   * ----------------------------------------------------------- */

  f->rho = c->rhou[dimension];

  if (c->rho > 0){
    float v = c->rhou[dimension] / c->rho;
    float p = GM1 * ( c->E - 0.5 * (c->rhou[0] * c->rhou[0] + c->rhou[1] * c->rhou[1] ) / c->rho );
    f->rhou[dimension] = c->rho * v * v + p;
    f->rhou[(dimension+1) % 2] = c->rhou[(dimension+1) % 2] * v;
    f->E = (c->E + p) * v;
  } else {
    f->rhou[0] = 0;
    f->rhou[1] = 0;
    f->E = 0;
  }
}




float gas_soundspeed(pstate* s){
  /*-----------------------------------------*/
  /* compute sound speed of ideal gas        */
  /*-----------------------------------------*/
  return sqrt(GAMMA * s->p / s->rho);
}




float gas_energy(pstate* s){
  /*-----------------------------------------*/
  /* compute total energy of a state         */
  /*-----------------------------------------*/
  return 0.5*s->rho * (s->u[0] * s->u[0] + s->u[1] * s->u[1]) + s->p/GM1;
}
