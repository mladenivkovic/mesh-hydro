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

#if LIMITER == NONE
  /* WAF calls limiter_get_phi directly, so you need to check here
   * whether that is the case. */
  limiter_get_phi_no_limiter(c, phi, dimension);
  return;
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
   * TODO: param documentation
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
   * TODO: param documentation
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
   * TODO: DOCUMENTATION
   * Compute the flow parameter r = ( u_i - u_i-1  ) / (u_i+1 - u_i) 
   * Also check whether you might compute junk by dividing by zero.
   * If you have, return something ridiculously high with the correct sign. 
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
  /* -------------------------------------------------
   * TODO: document!!!
   * ------------------------------------------------ */


  if (bottomleft == topleft){
    return ((topleft - topright)*1e18);
  } else{
    return ((topleft - topright)/(bottomleft - topleft));
  }
}
