/* Constant source acceleration         */

/* Written by Mladen Ivkovic, APR 2020
 * mladen.ivkovic@hotmail.com           */

#include "cell.h"
#include "params.h"

extern params pars;

void sources_get_acceleration(void) {
  /* --------------------------------------------
   * compute the accelleration for all cells
   * at the current time
   * -------------------------------------------- */

  /* If the acceleration is constant, don't compute it again. */
  if (pars.constant_acceleration_computed) { return; }

#if NDIM == 1
  for (int i = 0; i < pars.nxtot; i++) {
    cell* c   = &grid[i];
    c->acc[0] = pars.src_const_acc_x;
    c->acc[1] = pars.src_const_acc_y;
  }
#elif NDIM == 2
  for (int i = 0; i < pars.nxtot; i++) {
    for (int j = 0; j < pars.nxtot; j++) {
      cell* c   = &grid[i][j];
      c->acc[0] = pars.src_const_acc_x;
      c->acc[1] = pars.src_const_acc_y;
    }
  }
#endif

  /* if we have constant acceleration, set the flag that it's been dealt with */
  if (pars.constant_acceleration) { pars.constant_acceleration_computed = 1; }
}
