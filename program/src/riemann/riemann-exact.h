/* all stuff concerning exact Riemann solver */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef RIEMANN_EXACT_H
#define RIEMANN_EXACT_H
#include "gas.h"

void riemann_compute_star_states(oneDpstate *left, oneDpstate *right, float *pstar, float *ustar);
void riemann_sample_solution(oneDpstate* left, oneDpstate* right, float pstar, float ustar, oneDpstate* sol, float xovert, float* wavevel);
// int check_vacuum(oneDpstate *left, oneDpstate *right);
float fp(float pstar, oneDpstate *s, float A, float B, float a);
float dfpdp(float pstar, oneDpstate *s, float A, float B, float a);

#endif
