/* top level file for the limiter 
 * will be extended by the specified limiter*/


/* Written by Mladen Ivkovic, MAR 2020
 * mladen.ivkovic@hotmail.com           */


#include "gas.h"
#include "cell.h"
#include "limiter.h"
#include "params.h"

#include <math.h>
#include <stdio.h>


#if NDIM == 1
extern cell *grid;
#elif NDIM == 2
extern cell **grid;
#endif


extern params pars;




void limiter_get_phi(cell* c, pstate* phi, int dimension){
  /* ------------------------------------------------------------------------
   * Compute the flux limiter function phi_{i+1/2}
   *
   * cell* c:     for which cell i to work for
   * pstate* phi: where the limiter will be stored
   * dimension:   for which dimension we're working
   * ------------------------------------------------------------------------ */

#ifdef WAF 
#if LIMITER == NONE
  /* if we utilize a WAF method and no limiter, the implemented centered 
   * slope will give wrong results. Instead, catch it here and just return 
   * phi = 1. Then psi = 1 - (1 - |c|)phi(r) = |c| and we indeed get the
   * original method back. */
  phi->rho  = 1.;
  phi->u[0] = 1.;
  phi->u[1] = 1.;
  phi->p    = 1.;
  return;
#endif
#endif

  pstate Uim1, Ui, Uip1, Uip2;  /*U_i-1, U_i, U_i+1, U_i+2 */ 
  gas_init_pstate(&Uim1);
  gas_init_pstate(&Ui);
  gas_init_pstate(&Uip1);
  gas_init_pstate(&Uip2);


  /* Find which cells will come into play; */
  int i, j;
  cell_get_ij(c, &i, &j);

  float vel = 0;
#if NDIM == 1
  vel = c->prim.u[0];
  Uip2 = grid[i+2].prim;
  Uip1 = grid[i+1].prim;
  Ui = grid[i].prim;
  Uim1 = grid[i-1].prim;
#elif NDIM == 2
  if (dimension == 0){
    Uip2 = grid[i+2][j].prim;
    Uip1 = grid[i+1][j].prim;
    Ui = grid[i][j].prim;
    Uim1 = grid[i-1][j].prim;
  } else if (dimension == 1){
    vel = c->prim.u[1];
    Uip2 = grid[i][j+2].prim;
    Uip1 = grid[i][j+1].prim;
    Ui = grid[i][j].prim;
    Uim1 = grid[i][j-1].prim;
  }
#endif

  pstate r;
  gas_init_pstate(&r);
  limiter_get_r(&Uip2, &Uip1, &Ui, &Uim1, &r, vel);

  phi->rho  = limiter_phi_of_r(r.rho);
  phi->u[0] = limiter_phi_of_r(r.u[0]);
  phi->u[1] = limiter_phi_of_r(r.u[1]);
  phi->p    = limiter_phi_of_r(r.p);
}






void limiter_get_slope_left(cell* c, pstate* slope, int dimension){
  /* ------------------------------------------------------------------------
   * Compute the left slope of given cell c, i.e. the slope for the flux
   * F_{i-1/2}.
   * Just figure out which cell is one to the left and call 
   * limiter_get_slope_right for it instead of given cell c.
   *
   * cell* c:       cell for which we are working
   * pstate* slope: where the computed slope will be stored for all primitive
   *                states individually
   * int dimension: for which dimension we are working. 0: x, 1: y.
   * ------------------------------------------------------------------------ */

  int i, j;
  cell_get_ij(c, &i, &j);

  cell* left_cell;
#if NDIM == 1
  left_cell = &grid[i-1];
#elif NDIM == 2
  if (dimension == 0){
    left_cell = &grid[i-1][j];
  } else{
    left_cell = &grid[i][j-1];
  }
#endif

  limiter_get_slope_right(left_cell, slope, dimension);
}







void limiter_get_slope_right(cell* c, pstate* slope, int dimension){
  /* ------------------------------------------------------------------------
   * Compute the left slope of given cell c, i.e. the slope for the flux
   * F_{i+1/2}.
   * Remember: slope_i = 1/dx * phi(r_{i+1/2}) * (U_{i+1} - U_{i})
   *
   * cell* c:       cell for which we are working
   * pstate* slope: where the computed slope will be stored for all primitive
   *                states individually
   * int dimension: for which dimension we are working. 0: x, 1: y.
   * ------------------------------------------------------------------------ */

  pstate Uip1, Ui;
  gas_init_pstate(&Ui);
  gas_init_pstate(&Uip1);


  /* Find which cells will come into play*/
  int i, j;
  cell_get_ij(c, &i, &j);

  float vel = 0;
#if NDIM == 1
  vel = c->prim.u[0];
  Uip1 = grid[i+1].prim;
  Ui = grid[i].prim;
#elif NDIM == 2
  if (dimension == 0){
    vel = c->prim.u[0];
    Uip1 = grid[i+1][j].prim;
    Ui = grid[i][j].prim;
  } else if (dimension == 1){
    vel = c->prim.u[1];
    Uip1 = grid[i][j+1].prim;
    Ui = grid[i][j].prim;
  }
#endif

  /* Get the function phi */
  pstate phi;
  gas_init_pstate(&phi);
  limiter_get_phi(c, &phi, dimension);

  /* Now finally compute the actual slope */
  slope->rho  = phi.rho  * (Uip1.rho  - Ui.rho)  / pars.dx;
  slope->u[0] = phi.u[0] * (Uip1.u[0] - Ui.u[0]) / pars.dx;
  slope->u[1] = phi.u[1] * (Uip1.u[1] - Ui.u[1]) / pars.dx;
  slope->p    = phi.p    * (Uip1.p    - Ui.p)    / pars.dx;
}







void limiter_get_r(pstate* Uip2, pstate* Uip1, pstate* Ui, pstate* Uim1, pstate* r, float vel){
  /*----------------------------------------------------------------------------------------------
   * if v > 0: 
   *    compute r = (u_{i} - u_{i-1}) / (u_{i+1} - u_{i}) 
   * else:
   *    compute r = (u_{i+1} - u_{i+2}) / (u_{i} - u_{i+1})
   *
   * pstate* Uip2:  U_{i+2}
   * pstate* Uip1:  U_{i+1}
   * pstate* Ui:    U_{i}
   * pstate* Uim1:  U_{i-1}
   * pstate* r:     r for every primitive state
   * float vel:     advection velocity
   * ---------------------------------------------------------------------------------------------*/

  if (vel >= 0){
    r->rho =  limiter_r(Ui->rho,   Uim1->rho,  Uip1->rho);
    r->u[0] = limiter_r(Ui->u[0],  Uim1->u[0], Uip1->u[0]);
    r->u[1] = limiter_r(Ui->u[1],  Uim1->u[1], Uip1->u[1]);
    r->p =    limiter_r(Ui->p,     Uim1->p,    Uip1->p);
  } else{
    r->rho =  limiter_r(Uip1->rho,  Uip2->rho,  Ui->rho);
    r->u[0] = limiter_r(Uip1->u[0], Uip2->u[0], Ui->u[0]);
    r->u[1] = limiter_r(Uip1->u[1], Uip2->u[1], Ui->u[1]);
    r->p =    limiter_r(Uip1->p,    Uip2->p,    Ui->p);
  }
}






float limiter_r(float topleft, float topright, float bottomleft){
  /* --------------------------------------------------------------
   * if v > 0: 
   *    compute r = (u_{i} - u_{i-1}) / (u_{i+1} - u_{i}) 
   * else:
   *    compute r = (u_{i+1} - u_{i+2}) / (u_{i} - u_{i+1})
   *
   * In the tex documents, r for v < 0 is given as
   *    r = (u_{i+2} - u_{i+1}) / (u_{i+1} - u_{i})
   * which can be transformed into the above expression by multiplying
   * the numerator and the denominator by -1.
   * So then we can write
   *
   *          top_left - top_right
   * r = ---------------------------------
   *          bottom_left - top_left
   *
   * regardless of what sign the velocity v has. We only need to
   * switch what topleft, topright, and bottomleft are, which is
   * done in the function that is calling this one.
   * -------------------------------------------------------------- */


  if (bottomleft == topleft){
    return ((topleft - topright)*1e18);
  } else{
    return ((topleft - topright)/(bottomleft - topleft));
  }
}
