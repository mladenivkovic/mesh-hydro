/* First order Godunov scheme
 * declaration of functions additional to those in solver.h */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/

#ifndef GODUNOV_H
#define GODUNOV_H

#include "cell.h"

void solver_init_step(void);
void solver_compute_fluxes(float *dt, int dimension);
void solver_compute_cell_pair_flux(cell *left, cell *right, const float *dt, int dim);

#endif
