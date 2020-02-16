/* all stuff concerning exact Riemann solver */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef RIEMANN_EXACT_H
#define RIEMANN_EXACT_H
#include "gas.h"


extern int check_vacuum(pstate *left, pstate *right);
extern void compute_fluxes();
extern void compute_intercell_states();
extern void compute_star_pstate(pstate *left, pstate *right, pstate* starL, pstate* starR);
extern double compute_pstar(pstate *left, pstate *right);
extern double fp(double pguess, pstate *s, double gamma, double A, double B, double a);
extern double dfpdp(double pguess, pstate *s, double gamma, double A, double B, double a);
extern double fpfull(double fpL, double fpR, double delta_u);
extern double dfpdpfull(double dfpdpL, double dfpdpR);
extern double rho_star(pstate *s, pstate *star);
extern void compute_riemann(pstate* left, pstate* right, pstate* starL, pstate* starR, pstate* intercell);
extern void compute_riemann_vacuum(pstate* left, pstate* right, pstate* intercell);
extern double rho_fanL(pstate* s);
extern double u_fanL(pstate* s);
extern double p_fanL(pstate* s);
extern double rho_fanR(pstate* s);
extern double u_fanR(pstate* s);
extern double p_fanR(pstate* s);
#endif
