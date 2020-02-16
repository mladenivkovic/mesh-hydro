/* monotonized central difference limiter  */


/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */



#ifndef MC_H
#define MC_H


void limiter_get_r(pstate* Uip2, pstate* Uip1, pstate* Ui, pstate* Uim1, pstate* r, float vel);
float limiter_mc(float r);


#endif
