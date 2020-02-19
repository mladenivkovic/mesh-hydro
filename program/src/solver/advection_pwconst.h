/* Piecewise constant advection scheme 
 * declaration of functions additional to those in hydro.h */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */
/* -------------------------------------*/



#ifndef ADVECTION_H
#define ADVECTION_H

#ifdef ADVECTION_KEEP_VELOCITY_CONSTANT
void solver_advection_check_global_velocity();
#endif

void solver_init_step();
void solver_get_dt(MYFLOAT* dt);
void solver_compute_fluxes(int dimension);
void solver_advance_step(MYFLOAT* dt);

void solver_compute_cell_pair_flux(cell* c, cell* uw, cell* dw, int dimension);
void solver_update_state(cell* c, MYFLOAT dtdx);


#endif
