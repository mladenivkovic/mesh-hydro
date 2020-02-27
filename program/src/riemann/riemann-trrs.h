/* all stuff concerning Two Rarefaction Riemann solver */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef RIEMANN_TRRS_H
#define RIEMANN_TRRS_H
#include "gas.h"

void riemann_compute_star_states(pstate *left, pstate *right, float *pstar, float *ustar, int dimension);
void riemann_sample_solution(pstate* left, pstate* right, float pstar, float ustar, pstate* sol, float xovert, float* wavevel, int dim);

#endif
