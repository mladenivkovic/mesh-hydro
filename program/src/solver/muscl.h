/* MUSCL-Hancock scheme */
/* declaration of functions additional to those in solver.h */

/* Written by Mladen Ivkovic, APR 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/



#ifndef MUSCL_H
#define MUSCL_H

void solver_init_step();
void solver_compute_fluxes(float* dt, int dimension);
void solver_advance_step(float* dt, int dimension);

void solver_prepare_flux_computation(cell* c, float dthalf, int dim);
void solver_compute_cell_pair_flux(cell* left, cell* right, float* dt, int dim);
void solver_update_state(cell* left, cell* right, float dtdx);

#endif
