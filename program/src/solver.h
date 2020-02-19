/* top level file for the hydro method solver */


/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef SOLVER_H
#define SOLVER_H

#include "defines.h"


#if SOLVER == GODUNOV_UPWIND
#include "solver/godunov_upwind.h"
#elif SOLVER == ADVECTION_PWCONST
#include "solver/advection_pwconst.h"
#elif SOLVER == ADVECTION_PWLIN
#include "solver/advection_pwlin.h"
#endif


void solver_step(MYFLOAT* t, MYFLOAT* dt, int step,  int* write_output);

#endif
