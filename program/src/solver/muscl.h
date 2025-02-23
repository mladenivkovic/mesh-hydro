/* MUSCL-Hancock scheme */
/* declaration of functions additional to those in solver.h */

/* Written by Mladen Ivkovic, APR 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/

#ifndef MUSCL_H
#define MUSCL_H

#include "cell.h"

void solver_init_step(void);
void solver_compute_fluxes(float* dt, int dimension);
void solver_prepare_flux_computation(cell* c, float dthalf, int dim);
void solver_compute_cell_pair_flux(
  cell* left, cell* right, const float* dt, int dim
);

#endif
