/* all stuff concerning Riemann problem */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */


#ifndef RIEMANN_HLLC_H
#define RIEMANN_HLLC_H
#include "gas.h"

extern void compute_riemann(pstate* left, pstate* right, pstate* starL, pstate* starR, pstate* intercell);
extern void compute_fluxes();
extern void compute_wave_speeds(cstate left, cstate right,  double* SL, double* SR);
extern void Fhllc(double Sstar, double Sk, pstate wk, cstate uk, cstate fk, cstate* flux);
#endif
