/* monotonized central difference limiter  */


/* Written by Mladen Ivkovic, FEB 2020
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




void limiter_get_slope_right(cell* c, pstate* slope, int dimension){
  /* ------------------------------------------------------------------------
   * Compute the left slope of given cell c, i.e. the slope for the flux
   * F_{i+1/2}.
   * Remember: slope_i = 1/dx * phi(r_{i+1/2}) * (U_{i+1} - U_{i})
   * ------------------------------------------------------------------------ */
  int i, j;
  MYFLOAT vel = 0;
  pstate r, phi;
  /* cell left, right; */
  pstate Uim1, Ui, Uip1, Uip2;  /*U_i-1, U_i-2, U_i, U_i+1 */ 


  /* Find which cells wil come into play; This makes a difference for left/right */

  cell_get_ij(c, &i, &j);

  if (dimension == 0){
    vel = c->prim.ux;
#if NDIM == 1
    Uip2 = grid[i+2].prim;
    Uip1 = grid[i+1].prim;
    Ui = grid[i].prim;
    Uim1 = grid[i-1].prim;
#elif NDIM == 2
    Uip2 = grid[i+2][j].prim;
    Uip1 = grid[i+1][j].prim;
    Ui = grid[i][j].prim;
    Uim1 = grid[i-1][j].prim;
#endif
  } else if (dimension == 1){
    vel = c->prim.uy;
#if NDIM == 2
    Uip2 = grid[i][j+2].prim;
    Uip1 = grid[i][j+1].prim;
    Ui = grid[i][j].prim;
    Uim1 = grid[i][j-1].prim;
#endif
  }

  limiter_get_r(&Uip2, &Uip1, &Ui, &Uim1, &r, vel);


  phi.rho = limiter_mc(r.rho);
  phi.ux = limiter_mc(r.ux);
  phi.uy = limiter_mc(r.uy);
  phi.p = limiter_mc(r.p);

  /* Which Ui you take to do the difference below is different for left/right slope*/
  slope->rho = phi.rho * (Uip1.rho - Ui.rho) / pars.dx;
  slope->ux  = phi.ux  * (Uip1.ux  - Ui.ux)  / pars.dx;
  slope->uy  = phi.uy  * (Uip1.uy  - Ui.uy)  / pars.dx;
  slope->p   = phi.p   * (Uip1.p   - Ui.p)   / pars.dx;
}





void limiter_get_slope_left(cell* c, pstate* slope, int dimension){
  /* ------------------------------------------------------------------------
   * Compute the left slope of given cell c, i.e. the slope for the flux
   * F_{i-1/2}.
   * Remember: slope_{i-1} = 1/dx * phi(r_{i-1/2}) * (U_i - U_{i-1})
   * ------------------------------------------------------------------------ */
  int i, j;
  MYFLOAT vel = 0;
  pstate r, phi;
  pstate Uim1, Uim2, Ui, Uip1;  /*U_i-1, U_i-2, U_i, U_i+1 */


  /* Find which cells wil come into play; This makes a difference for left/right */

  cell_get_ij(c, &i, &j);

  if (dimension == 0){
    vel = c->prim.ux;
#if NDIM == 1
    Uip1 = grid[i+1].prim;
    Ui = grid[i].prim;
    Uim1 = grid[i-1].prim;
    Uim2 = grid[i-2].prim;
#elif NDIM == 2
    Uip1 = grid[i+1][j].prim;
    Ui = grid[i][j].prim;
    Uim1 = grid[i-1][j].prim;
    Uim2 = grid[i-2][j].prim;
#endif
  } else if (dimension == 1){
    vel = c->prim.uy;
#if NDIM == 2
    Uip1 = grid[i][j+1].prim;
    Ui = grid[i][j].prim;
    Uim1 = grid[i][j-1].prim;
    Uim2 = grid[i][j-2].prim;
#endif
  }

  limiter_get_r(&Uip1, &Ui, &Uim1, &Uim2, &r, vel);

  phi.rho = limiter_mc(r.rho);
  phi.ux = limiter_mc(r.ux);
  phi.uy = limiter_mc(r.uy);
  phi.p = limiter_mc(r.p);

  /* Which Ui you take to do the difference below is different for left/right slope*/
  slope->rho = phi.rho * (Ui.rho - Uim1.rho) / pars.dx;
  slope->ux  = phi.ux  * (Ui.ux  - Uim1.ux)  / pars.dx;
  slope->uy  = phi.uy  * (Ui.uy  - Uim1.uy)  / pars.dx;
  slope->p   = phi.p   * (Ui.p   - Uim1.p)   / pars.dx;
}








void limiter_get_r(pstate* Uip1, pstate* Ui, pstate* Uim1, pstate* Uim2, pstate* r, MYFLOAT vel){

  /* Also check whether you might compute junk by dividing by zero.
   * If you have, return something ridiculously high with the correct sign. */

  if (vel >= 0){
    if (Ui->rho == Uim1->rho){
      r->rho = (Uim1->rho - Uim2->rho)*1e18;
    } else {
      r->rho = (Uim1->rho - Uim2->rho)/(Ui->rho - Uim1->rho);
    }
    if (Ui->ux == Uim1->ux){
      r->ux = (Uim1->ux - Uim2->ux)*1e18;
    } else {
      r->ux = (Uim1->ux - Uim2->ux)/(Ui->ux - Uim1->ux);
    }
    if (Ui->uy == Uim1->uy){
      r->uy = (Uim1->uy - Uim2->uy)*1e18;
    } else {
      r->uy = (Uim1->uy - Uim2->uy)/(Ui->uy - Uim1->uy);
    }
    if (Ui->p == Uim1->p){
      r->p = (Uim1->p - Uim2->p)*1e18;
    } else {
      r->p = (Uim1->p - Uim2->p)/(Ui->p - Uim1->p);
    }
  } else { /* vel < 0 */
    if (Ui->rho == Uim1->rho){
      r->rho = (Uip1->rho - Ui->rho)*1e18;
    } else {
      r->rho = (Uip1->rho - Ui->rho)/(Ui->rho - Uim1->rho);
    }
    if (Ui->ux == Uim1->ux){
      r->ux = (Uip1->ux - Ui->ux)*1e18;
    } else{
      r->ux = (Uip1->ux - Ui->ux)/(Ui->ux - Uim1->ux);
    }
    if (Ui->uy == Uim1->uy){
      r->uy = (Uip1->uy - Ui->uy)*1e18;
    } else{
      r->uy = (Uip1->uy - Ui->uy)/(Ui->uy - Uim1->uy);
    }
    if (Ui->p == Uim1->p){
      r->p = (Uip1->p - Ui->p)*1e18;
    } else {
      r->p = (Uip1->p - Ui->p)/(Ui->p - Uim1->p);
    }
  }
}






MYFLOAT limiter_mc(MYFLOAT r){
  /* -----------------------------------------------------
   * Computes the phi(r) for the superbee slope limiter
   * ----------------------------------------------------- */

  MYFLOAT min = 0.0;
  MYFLOAT max = 0.0;

  /* min((1+r)/2, 2) */
  MYFLOAT temp = 0.5*(1+r);
  if (temp > 2.0) temp = 2.0;
  min = temp;

  /* min( min((1+r)/2, 2), 2r) */
  temp = 2 * r;
  if (temp > min) temp = min;       
  min = temp;
  
  /* max( 0, min(...) ) */
  if (min > max) max = min;

  return(max);
}
