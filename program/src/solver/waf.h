/* First order Godunov scheme */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/



#ifndef WAF_H
#define WAF_H

void solver_init_step();
void solver_get_dt(float* dt, int step);
void solver_compute_fluxes(float* dt, int dimension);
void solver_advance_step(float* dt, int dimension);

void solver_prepare_flux_computation(cell* left, cell* right, int dim);
void solver_compute_cell_pair_flux(cell* left, cell* right, float* dt, int dim);
void solver_update_state(cell* left, cell* right, float dtdx);

#endif