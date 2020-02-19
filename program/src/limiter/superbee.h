/* SUPERBEE limiter  */


/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */



#ifndef SUPERBEE_H
#define SUPERBEE_H


void limiter_get_r(pstate* Uip2, pstate* Uip1, pstate* Ui, pstate* Uim1, pstate* r, MYFLOAT vel);
MYFLOAT limiter_superbee(MYFLOAT r);


#endif
