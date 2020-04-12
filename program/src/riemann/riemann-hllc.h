/* all stuff concerning Riemann problem */
/* This file extends and is included by /program/src/riemann.h */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */


#ifndef RIEMANN_HLLC_H
#define RIEMANN_HLLC_H
#include "gas.h"

void riemann_solve_hllc(pstate* left, pstate* right, cstate* sol, 
    float xovert, int dimension);

void riemann_compute_wave_speed_estimates(pstate* left, pstate* right, 
    float* SL, float* SR, float* Sstar, int dim);

float qLR(float pstar, float pLR);

void riemann_sample_hllc_solution(pstate* left, pstate* right, float SL, 
    float SR, float Sstar, cstate* sol, float xovert, int dim);

void riemann_get_hllc_full_solution(pstate* left, pstate* right, 
    float* S, cstate* fluxes, float* delta_q, int dim);

#endif
