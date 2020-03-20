/* Piecewise linear advection scheme 
 * declaration of functions additional to those in solver.h */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/



#ifndef ADVECTION_H
#define ADVECTION_H

void solver_init_step();
void solver_compute_fluxes(float* dt, int dimension);
void solver_advance_step(float* dt);

void solver_compute_cell_pair_flux(cell* c, cell* uw, cell* dw, float* dt, int dim);
void solver_compute_slope(cell* uw, cell* dw,  pstate* slope);
#endif
