/* Weighted Average Flux ADVECTION scheme
 * declaration of functions additional to those in solver.h */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/

#ifndef ADVECTION_WAF_H
#define ADVECTION_WAF_H

#include "cell.h"

#ifdef ADVECTION_KEEP_VELOCITY_CONSTANT
void solver_advection_check_global_velocity(void);
#endif

void solver_init_step(void);
void solver_compute_fluxes(float* dt, int dimension);
void solver_compute_cell_pair_flux(cell* c, cell* n, const float* dt, int dim);

#endif
