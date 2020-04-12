/* top level file for the limiter */


/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef LIMITER_H
#define LIMITER_H

#include "defines.h"

#if LIMITER == NONE
#include "limiter/no_limiter.h"
#elif LIMITER == MINMOD
#include "limiter/minmod.h"
#elif LIMITER == SUPERBEE
#include "limiter/superbee.h"
#elif LIMITER == VANLEER
#include "limiter/van_leer.h"
#elif LIMITER == MC
#include "limiter/monotonized_central_difference.h"
#endif


void limiter_get_slope_left(cell* c, pstate* slope, int dimension);
void limiter_get_slope_right(cell* c, pstate* slope, int dimension);
void limiter_get_phi(cell* c, pstate* phi, int dimension);
void limiter_get_r(pstate* Uip2, pstate* Uip1, pstate* Ui, pstate* Uim1, pstate* r, float vel);
float limiter_phi_of_r(float r);
float limiter_r(float topleft, float topright, float bottomleft);

#endif
