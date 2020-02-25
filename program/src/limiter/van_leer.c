/* VAN_LEER limiter  */


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
  float vel = 0;
  pstate r, phi;
  gas_init_pstate(&r);
  gas_init_pstate(&phi);

  pstate Uim1, Ui, Uip1, Uip2;  /*U_i-1, U_i-2, U_i, U_i+1 */ 
  gas_init_pstate(&Uim1);
  gas_init_pstate(&Ui);
  gas_init_pstate(&Uip1);
  gas_init_pstate(&Uip2);


  /* Find which cells wil come into play; This makes a difference for left/right */

  cell_get_ij(c, &i, &j);

  if (dimension == 0){
    vel = c->prim.u[0];
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
    vel = c->prim.u[1];
#if NDIM == 2
    Uip2 = grid[i][j+2].prim;
    Uip1 = grid[i][j+1].prim;
    Ui = grid[i][j].prim;
    Uim1 = grid[i][j-1].prim;
#endif
  }

  limiter_get_r(&Uip2, &Uip1, &Ui, &Uim1, &r, vel);


  phi.rho = limiter_vanleer(r.rho);
  phi.u[0] = limiter_vanleer(r.u[0]);
  phi.u[1] = limiter_vanleer(r.u[1]);
  phi.p = limiter_vanleer(r.p);

  /* Which Ui you take to do the difference below is different for left/right slope*/
  slope->rho = phi.rho * (Uip1.rho - Ui.rho) / pars.dx;
  slope->u[0]  = phi.u[0]  * (Uip1.u[0]  - Ui.u[0])  / pars.dx;
  slope->u[1]  = phi.u[1]  * (Uip1.u[1]  - Ui.u[1])  / pars.dx;
  slope->p   = phi.p   * (Uip1.p   - Ui.p)   / pars.dx;
}





void limiter_get_slope_left(cell* c, pstate* slope, int dimension){
  /* ------------------------------------------------------------------------
   * Compute the left slope of given cell c, i.e. the slope for the flux
   * F_{i-1/2}.
   * Remember: slope_{i-1} = 1/dx * phi(r_{i-1/2}) * (U_i - U_{i-1})
   * ------------------------------------------------------------------------ */
  int i, j;
  float vel = 0;
  pstate r, phi;
  gas_init_pstate(&r);
  gas_init_pstate(&phi);

  pstate Uim1, Uim2, Ui, Uip1;  /*U_i-1, U_i-2, U_i, U_i+1 */
  gas_init_pstate(&Uim1);
  gas_init_pstate(&Uim2);
  gas_init_pstate(&Ui);
  gas_init_pstate(&Uip1);


  /* Find which cells wil come into play; This makes a difference for left/right */

  cell_get_ij(c, &i, &j);

  if (dimension == 0){
    vel = c->prim.u[0];
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
    vel = c->prim.u[1];
#if NDIM == 2
    Uip1 = grid[i][j+1].prim;
    Ui = grid[i][j].prim;
    Uim1 = grid[i][j-1].prim;
    Uim2 = grid[i][j-2].prim;
#endif
  }

  limiter_get_r(&Uip1, &Ui, &Uim1, &Uim2, &r, vel);

  phi.rho = limiter_vanleer(r.rho);
  phi.u[0] = limiter_vanleer(r.u[0]);
  phi.u[1] = limiter_vanleer(r.u[1]);
  phi.p = limiter_vanleer(r.p);

  /* Which Ui you take to do the difference below is different for left/right slope*/
  slope->rho = phi.rho * (Ui.rho - Uim1.rho) / pars.dx;
  slope->u[0]  = phi.u[0]  * (Ui.u[0]  - Uim1.u[0])  / pars.dx;
  slope->u[1]  = phi.u[1]  * (Ui.u[1]  - Uim1.u[1])  / pars.dx;
  slope->p   = phi.p   * (Ui.p   - Uim1.p)   / pars.dx;
}






void limiter_get_r(pstate* Uip1, pstate* Ui, pstate* Uim1, pstate* Uim2, pstate* r, float vel){

  /* Also check whether you might compute junk by dividing by zero.
   * If you have, return something ridiculously high with the correct sign. */

  if (vel >= 0){
    if (Ui->rho == Uim1->rho){
      r->rho = (Uim1->rho - Uim2->rho)*1e18;
    } else {
      r->rho = (Uim1->rho - Uim2->rho)/(Ui->rho - Uim1->rho);
    }
    if (Ui->u[0] == Uim1->u[0]){
      r->u[0] = (Uim1->u[0] - Uim2->u[0])*1e18;
    } else {
      r->u[0] = (Uim1->u[0] - Uim2->u[0])/(Ui->u[0] - Uim1->u[0]);
    }
    if (Ui->u[1] == Uim1->u[1]){
      r->u[1] = (Uim1->u[1] - Uim2->u[1])*1e18;
    } else {
      r->u[1] = (Uim1->u[1] - Uim2->u[1])/(Ui->u[1] - Uim1->u[1]);
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
    if (Ui->u[0] == Uim1->u[0]){
      r->u[0] = (Uip1->u[0] - Ui->u[0])*1e18;
    } else{
      r->u[0] = (Uip1->u[0] - Ui->u[0])/(Ui->u[0] - Uim1->u[0]);
    }
    if (Ui->u[1] == Uim1->u[1]){
      r->u[1] = (Uip1->u[1] - Ui->u[1])*1e18;
    } else{
      r->u[1] = (Uip1->u[1] - Ui->u[1])/(Ui->u[1] - Uim1->u[1]);
    }
    if (Ui->p == Uim1->p){
      r->p = (Uip1->p - Ui->p)*1e18;
    } else {
      r->p = (Uip1->p - Ui->p)/(Ui->p - Uim1->p);
    }
  }
}






float limiter_vanleer(float r){
  /* -----------------------------------------------------
   * Computes the phi(r) for the superbee slope limiter
   * ----------------------------------------------------- */

  float phi = (r + fabs(r))/(1 + fabs(r));
  return(phi);
}
