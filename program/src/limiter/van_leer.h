/* VAN_LEER limiter  */


/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */



#ifndef VAN_LEER_H
#define VAN_LEER_H


void limiter_get_r(pstate* Uip2, pstate* Uip1, pstate* Ui, pstate* Uim1, pstate* r, MYFLOAT vel);
MYFLOAT limiter_vanleer(MYFLOAT r);


#endif
