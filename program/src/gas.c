/* ideal gas related stuff */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#include <stdio.h>
#include <math.h>

#include "gas.h"
#include "defines.h"



void gas_init_pstate(pstate* s){
  /*-------------------------------------------------*/
  /* This function sets the pstate to zero           */
  /*-------------------------------------------------*/

  s->rho  = 0;
  s->ux   = 0;
  s->uy   = 0;
  s->p    = 0;

} 

void gas_init_cstate(cstate *s){
  /*-------------------------------------------------*/
  /* This function sets the cstate to zero           */
  /*-------------------------------------------------*/

  s->rho   = 0;
  s->rhoux = 0;
  s->rhouy = 0;
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
  return 0.5*s->rho * (s->ux * s->ux + s->uy * s->uy) + s->p/GM1;
}
