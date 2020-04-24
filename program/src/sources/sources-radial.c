/* radial constant source acceleration         */
/* Definition: Positive acceleration points from
 * the center of the box (or line in 1D) outwards */

/* Written by Mladen Ivkovic, APR 2020
 * mladen.ivkovic@hotmail.com           */

#include "cell.h"
#include "defines.h"
#include "params.h"

#include <math.h>
#include <stdio.h>

#if NDIM == 1
extern cell *grid;
#elif NDIM == 2
extern cell **grid;
#endif

extern params pars;


void sources_get_acceleration(){
  /* --------------------------------------------
   * compute the accelleration for all cells
   * at the current time
   * -------------------------------------------- */

  /* If the acceleration is constant, don't compute it again. */
  if (pars.constant_acceleration_computed) return;



#if NDIM == 1

  float cx = 0.5 * BOXLEN; /* x coordinate of center of box */

  for (int i = 0; i < pars.nxtot; i++){
    cell* c = &grid[i];

    float s = 1;
    if (c->x < cx) s = -1;

    c->acc[0] = s * pars.src_const_acc_r;
    c->acc[1] = 0.;
  }


#elif NDIM == 2


  float cx = 0.5 * BOXLEN; /* x coordinate of center of box */
  float cy = 0.5 * BOXLEN; /* y coordinate of center of box */

  for (int j = 0; j < pars.nxtot; j++){
    for (int i = 0; i < pars.nxtot; i++){
      cell* c = &grid[i][j];

      /* pretend center of box is origin */
      float x = c->x;
      float y = c->y;
      float dx = x - cx;
      float dy = y - cy;
      float alpha = atan(dy/dx);
      float s = 1;
      if (dx < 0) s = -1;

      c->acc[0] = cos(alpha) * pars.src_const_acc_r * s;
      c->acc[1] = sin(alpha) * pars.src_const_acc_r * s;
    }
  }
#endif

  /* if we have constant acceleration, set the flag that it's been dealt with */
  if (pars.constant_acceleration) pars.constant_acceleration_computed = 1;
}
