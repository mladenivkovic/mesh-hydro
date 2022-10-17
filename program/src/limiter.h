/* top level file for the limiter */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef LIMITER_H
#define LIMITER_H

#include "defines.h"
#include "cell.h"

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

void limiter_get_limited_slope(cell *c, cstate *slope, int dimension);
void limiter_get_advection_slope_left(cell *c, pstate *slope, int dimension);
void limiter_get_advection_slope_right(cell *c, pstate *slope, int dimension);
void limiter_get_phi(cell *, pstate *phi, int dimension);
void limiter_get_r_pstate(pstate *Uip2, pstate *Uip1, pstate *Ui, pstate *Uim1,
                          pstate *r, float vel);
void limiter_get_r_cstate(cstate *Uip1, cstate *Ui, cstate *Uim1, cstate *r);
float limiter_r(float topleft, float topright, float bottomleft);

/* declared here, should be defined in individual files */
float limiter_phi_of_r(float r);
float limiter_xi_of_r(float r);

#endif
