/* Global source term related functions */

/* Written by Mladen Ivkovic, APR 2020
 * mladen.ivkovic@hotmail.com           */

#include "sources.h"

#include "cell.h"
#include "gas.h"
#include "integrate.h"
#include "params.h"

extern params pars;

/**
 * Compute the source vector S for a given cell we assume that the CONSERVED
 * quantities are up to date
 * TODO: dox
 */
void sources_get_source_vector(cstate* s, const float acc[2], cstate* cons) {

  s->rho     = 0.;
  s->rhou[0] = cons->rho * acc[0];
  s->rhou[1] = cons->rho * acc[1];
  s->E       = cons->rhou[0] * acc[0] + cons->rhou[1] * acc[1];
}


/**
 * Compute the sources update over given time step dt for all cells
 */
void sources_update_state(float dt) {

  /* first, compute the current source terms */
  sources_get_acceleration();

#if NDIM == 1
  for (int i = 0; i < pars.nxtot; i++) {
    cell* c = &grid[i];
    integrate(&c->cons, c->acc, dt, &c->cons);
  }
#elif NDIM == 2
  for (int i = 0; i < pars.nxtot; i++) {
    for (int j = 0; j < pars.nxtot; j++) {
      cell* c = &grid[i][j];
      integrate(&c->cons, c->acc, dt, &c->cons);
    }
  }
#endif
}
