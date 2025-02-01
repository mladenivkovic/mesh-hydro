/* radial constant source acceleration         */
/* Definition: Positive acceleration points from
 * the center of the box (or line in 1D) outwards */

/* Written by Mladen Ivkovic, APR 2020
 * mladen.ivkovic@hotmail.com           */

#include <math.h>

#include "cell.h"
#include "params.h"

extern params pars;

void sources_get_acceleration(void) {
  /* --------------------------------------------
   * compute the accelleration for all cells
   * at the current time
   * -------------------------------------------- */

  /* If the acceleration is constant, don't compute it again. */
  if (pars.constant_acceleration_computed) {
    return;
  }

#if NDIM == 1

  float cx = 0.5 * BOXLEN; /* x coordinate of center of box */

  for (int i = 0; i < pars.nxtot; i++) {
    cell *c = &grid[i];

    float s = 1.f;
    if (c->x < cx) {
      s = -1.f;
    }

    c->acc[0] = s * pars.src_const_acc_r;
    c->acc[1] = 0.f;
  }

#elif NDIM == 2

  float cx = 0.5f * BOXLEN; /* x coordinate of center of box */
  float cy = 0.5f * BOXLEN; /* y coordinate of center of box */

  for (int j = 0; j < pars.nxtot; j++) {
    for (int i = 0; i < pars.nxtot; i++) {
      cell *c = &grid[i][j];

      /* pretend center of box is origin */
      float x = c->x;
      float y = c->y;
      float dx = x - cx;
      float dy = y - cy;
      float alpha = atanf(dy / dx);
      float s = 1;
      if (dx < 0) {
        s = -1.f;
      }

      c->acc[0] = cosf(alpha) * pars.src_const_acc_r * s;
      c->acc[1] = sinf(alpha) * pars.src_const_acc_r * s;
    }
  }
#endif

  /* if we have constant acceleration, set the flag that it's been dealt with */
  if (pars.constant_acceleration) {
    pars.constant_acceleration_computed = 1;
  }
}
