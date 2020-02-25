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
  s->u = malloc(2 * sizeof(float));
  s->u[0] = 0;
  s->u[1] = 0;
  s->p    = 0;

} 

void gas_init_cstate(cstate *s){
  /*-------------------------------------------------*/
  /* This function sets the cstate to zero           */
  /*-------------------------------------------------*/

  s->rho   = 0;
  s->rhou = malloc(2 * sizeof(float));
  s->rhou[0] = 0;
  s->rhou[1] = 0;
  s->E     = 0;

}





float gas_soundspeed(float p, float rho){
  /*-----------------------------------------*/
  /* compute sound speed of ideal gas        */
  /*-----------------------------------------*/
  return sqrtf(GAMMA * p / rho);
}




float gas_energy(pstate* s){
  /*-----------------------------------------*/
  /* compute total energy of a state         */
  /*-----------------------------------------*/
  return 0.5*s->rho * (s->u[0] * s->u[0] + s->u[1] * s->u[1]) + s->p/GM1;
}
