/* top level file for the Riemann solver */


/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef RIEMANN_H
#define RIEMANN_H

#include "defines.h"

#if RIEMANN == EXACT
#include "riemann/riemann-exact.h"
#elif RIEMANN == TRRS
#include "riemann/riemann-trrs.h"
#elif RIEMANN == TSRS
#include "riemann/riemann-tsrs.h"
#elif RIEMANN == HLLC
#include "riemann/riemann-hllc.h"
#endif

void riemann_solve(pstate* left, pstate* right, pstate* sol, 
    float xovert, int dimension);

int riemann_has_vacuum(pstate *left, pstate *right, int dimension);

void riemann_compute_vacuum_solution(pstate* left, pstate* right, 
    pstate* sol, float xovert, int dim);

void riemann_compute_star_states(pstate *left, pstate *right, 
    float *pstar, float *ustar, int dimension);

void riemann_sample_solution(pstate* left, pstate* right, 
    float pstar, float ustar, pstate* sol, float xovert, int dim);

void riemann_get_full_solution(pstate* left, pstate* right, 
    float S[3], cstate fluxes[4], float delta_q[3], int dim);

void riemann_get_full_vacuum_solution(pstate* left, pstate* right, 
    float S[3], cstate fluxes[4], float delta_q[3], int dim);
#endif
