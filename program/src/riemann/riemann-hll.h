/* all stuff concerning Riemann problem */


/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef RIEMANN_HLL_H
#define RIEMANN_HLL_H
#include "gas.h"

extern void compute_riemann(pstate* left, pstate* right, pstate* starL, pstate* starR, pstate* intercell);
extern void compute_fluxes();
// extern void compute_uhll(cstate* left, cstate* right, cstate* u_hll, double SL, double SR);
extern void compute_wave_speeds(cstate left, cstate right,  double* SL, double* SR);
extern double Fhll(double ul, double ur, double fl, double fr, double SL, double SR);
#endif
