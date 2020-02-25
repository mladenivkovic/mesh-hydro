/* all stuff concerning exact Riemann solver */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef RIEMANN_EXACT_H
#define RIEMANN_EXACT_H
#include "gas.h"

void riemann_compute_star_states(pstate *left, pstate *right, float *pstar, float *ustar, int dimension);
void riemann_sample_solution(pstate* left, pstate* right, float pstar, float ustar, pstate* sol, float xovert, float* wavevel, int dim);
// int check_vacuum(pstate *left, pstate *right);
float fp(float pstar, pstate *s, float A, float B, float a);
float dfpdp(float pstar, pstate *s, float A, float B, float a);

#endif
