/* top level file for the hydro method solver */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef SOLVER_H
#define SOLVER_H

#include "cell.h"
#include "defines.h"

#if SOLVER == GODUNOV
#include "solver/godunov.h"
#elif SOLVER == WAF
#include "solver/waf.h"
#elif SOLVER == MUSCL
#include "solver/muscl.h"
#elif SOLVER == ADVECTION_PWCONST
#include "solver/advection_pwconst.h"
#elif SOLVER == ADVECTION_PWLIN
#include "solver/advection_pwlin.h"
#elif SOLVER == ADVECTION_WAF
#include "solver/advection_waf.h"
#endif

void solver_step(const float *t, float *dt, int step, int *write_output);

void solver_get_advection_dt(float *dt);
void solver_advance_step_advection(const float *dt);
void solver_update_state_advection(cell *c, float dtdx);

void solver_get_hydro_dt(float *dt, int step);
void solver_advance_step_hydro(const float *dt, int dimension);
void solver_update_state_hydro(cell *left, cell *right, float dtdx);

#ifdef ADVECTION_KEEP_VELOCITY_CONSTANT
void solver_advection_check_global_velocity(void);
#endif

#endif
