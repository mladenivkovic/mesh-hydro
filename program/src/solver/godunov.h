/* First order Godunov scheme */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/



#ifndef GODUNOV_H
#define GODUNOV_H

void solver_init_step();
void solver_get_dt(float* dt);
void solver_compute_fluxes(float* dt, int dimension);
void solver_advance_step(float* dt, int dimension);

void solver_compute_cell_pair_flux(cell* left, cell* right, float* dt, int dim);
void solver_update_state(cell* left, cell* right, float dtdx);

#endif
