/* top level file for the limiter
 * will be extended by the specified limiter*/

/* Written by Mladen Ivkovic, MAR 2020
 * mladen.ivkovic@hotmail.com           */

#include "limiter.h"
#include "cell.h"
#include "gas.h"
#include "params.h"

#include <math.h>
#include <stdio.h>

extern params pars;

/**
 * Compute the slope of given cell c using slope limiters, as it is needed
 * for the MUSCL-Hancock scheme.
 * Remember:
 *  slope = 0.5(1 + omega)(U_{i} - U_{i-1}) + 0.5(1 - omega)(U_{i+1} - U_{i})
 * where omega is set in defines.h
 *
 * @param cell* c:       cell for which we are working
 * @param cstate* slope: where the computed slope will be stored for all
 *                       conservative states individually
 * @param int dimension: for which dimension we are working. 0: x, 1: y.
 */
void limiter_get_limited_slope(cell *c, cstate *slope, int dimension) {

  cstate Uip1, Ui, Uim1;
  gas_init_cstate(&Ui);
  gas_init_cstate(&Uip1);
  gas_init_cstate(&Uim1);

  /* Find which cells will come into play*/
  int i, j;
  cell_get_ij(c, &i, &j);

#if NDIM == 1
  Uip1 = grid[i + 1].cons;
  Ui = grid[i].cons;
  Uim1 = grid[i - 1].cons;
#elif NDIM == 2
  if (dimension == 0) {
    Uip1 = grid[i + 1][j].cons;
    Ui = grid[i][j].cons;
    Uim1 = grid[i - 1][j].cons;
  } else if (dimension == 1) {
    Uip1 = grid[i][j + 1].cons;
    Ui = grid[i][j].cons;
    Uim1 = grid[i][j - 1].cons;
  }
#endif

  /* Get the slope limiter xi */
  cstate r;
  gas_init_cstate(&r);
  limiter_get_r_cstate(&Uip1, &Ui, &Uim1, &r);

  cstate xi;
  gas_init_cstate(&xi);

  xi.rho = limiter_xi_of_r(r.rho);
  xi.rhou[0] = limiter_xi_of_r(r.rhou[0]);
  xi.rhou[1] = limiter_xi_of_r(r.rhou[1]);
  xi.E = limiter_xi_of_r(r.E);

  /* Now finally compute the actual slope */
  slope->rho =
      xi.rho * 0.5 *
      ((1. + OMEGA) * (Ui.rho - Uim1.rho) + (1. - OMEGA) * (Uip1.rho - Ui.rho));
  slope->rhou[0] = xi.rhou[0] * 0.5 *
                   ((1. + OMEGA) * (Ui.rhou[0] - Uim1.rhou[0]) +
                    (1. - OMEGA) * (Uip1.rhou[0] - Ui.rhou[0]));
  slope->rhou[1] = xi.rhou[1] * 0.5 *
                   ((1. + OMEGA) * (Ui.rhou[1] - Uim1.rhou[1]) +
                    (1. - OMEGA) * (Uip1.rhou[1] - Ui.rhou[1]));
  slope->E = xi.E * 0.5 *
             ((1. + OMEGA) * (Ui.E - Uim1.E) + (1. - OMEGA) * (Uip1.E - Ui.E));
}


/**
 * Compute the flux limiter function phi_{i+1/2}
 *
 * @param cell* c:     for which cell i to work for
 * @param pstate* phi: where the limiter will be stored
 * @param dimension:   for which dimension we're working
 */
void limiter_get_phi(cell *c, pstate *phi, int dimension) {

#if (SOLVER == ADVECTION_WAF) && (LIMITER == NONE)
  /* if we utilize a WAF method and no limiter, the implemented centered
   * slope will give wrong results. Instead, catch it here and just return
   * phi = 1. Then psi = 1 - (1 - |c|)phi(r) = |c| and we indeed get the
   * original method back. */
  phi->rho = 1.;
  phi->u[0] = 1.;
  phi->u[1] = 1.;
  phi->p = 1.;
  return;
#endif

  pstate Uim1, Ui, Uip1, Uip2; /*U_i-1, U_i, U_i+1, U_i+2 */
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
  Uip2 = grid[i + 2].prim;
  Uip1 = grid[i + 1].prim;
  Ui = grid[i].prim;
  Uim1 = grid[i - 1].prim;
#elif NDIM == 2
  vel = c->prim.u[dimension];
  if (dimension == 0) {
    Uip2 = grid[i + 2][j].prim;
    Uip1 = grid[i + 1][j].prim;
    Ui = grid[i][j].prim;
    Uim1 = grid[i - 1][j].prim;
  } else if (dimension == 1) {
    Uip2 = grid[i][j + 2].prim;
    Uip1 = grid[i][j + 1].prim;
    Ui = grid[i][j].prim;
    Uim1 = grid[i][j - 1].prim;
  }
#endif

  pstate r;
  gas_init_pstate(&r);
  limiter_get_r_pstate(&Uip2, &Uip1, &Ui, &Uim1, &r, vel);

  phi->rho = limiter_phi_of_r(r.rho);
  phi->u[0] = limiter_phi_of_r(r.u[0]);
  phi->u[1] = limiter_phi_of_r(r.u[1]);
  phi->p = limiter_phi_of_r(r.p);
}


/**
 * Compute the left slope of given cell c, i.e. the slope for the flux
 * F_{i-1/2}.
 * Just figure out which cell is one to the left and call
 * limiter_get_slope_right for it instead of given cell c.
 *
 * @param cell* c:       cell for which we are working
 * @param pstate* slope: where the computed slope will be stored for all
 *                       primitive states individually
 * @param int dimension: for which dimension we are working. 0: x, 1: y.
 */
void limiter_get_advection_slope_left(cell *c, pstate *slope, int dimension) {

  int i, j;
  cell_get_ij(c, &i, &j);

  cell *left_cell;
#if NDIM == 1
  left_cell = &grid[i - 1];
#elif NDIM == 2
  if (dimension == 0) {
    left_cell = &grid[i - 1][j];
  } else {
    left_cell = &grid[i][j - 1];
  }
#endif

  limiter_get_advection_slope_right(left_cell, slope, dimension);
}


/**
 * Compute the right slope of given cell c, i.e. the slope for the flux
 * F_{i+1/2}.
 * Remember: slope_i = 1/dx * phi(r_{i+1/2}) * (U_{i+1} - U_{i})
 *
 * @param cell* c:       cell for which we are working
 * @param pstate* slope: where the computed slope will be stored for all
 *                       primitive states individually
 * @param int dimension: for which dimension we are working. 0: x, 1: y.
 */
void limiter_get_advection_slope_right(cell *c, pstate *slope, int dimension) {

  pstate Uip1, Ui;
  gas_init_pstate(&Ui);
  gas_init_pstate(&Uip1);

  /* Find which cells will come into play*/
  int i, j;
  cell_get_ij(c, &i, &j);

#if NDIM == 1
  Uip1 = grid[i + 1].prim;
  Ui = grid[i].prim;
#elif NDIM == 2
  if (dimension == 0) {
    Uip1 = grid[i + 1][j].prim;
    Ui = grid[i][j].prim;
  } else if (dimension == 1) {
    Uip1 = grid[i][j + 1].prim;
    Ui = grid[i][j].prim;
  }
#endif

  /* Get the function phi */
  pstate phi;
  gas_init_pstate(&phi);
  limiter_get_phi(c, &phi, dimension);

  /* Now finally compute the actual slope */
  slope->rho = phi.rho * (Uip1.rho - Ui.rho) / pars.dx;
  slope->u[0] = phi.u[0] * (Uip1.u[0] - Ui.u[0]) / pars.dx;
  slope->u[1] = phi.u[1] * (Uip1.u[1] - Ui.u[1]) / pars.dx;
  slope->p = phi.p * (Uip1.p - Ui.p) / pars.dx;
}


/**
 * Compute the flow parameter r for every component of the primitive states.
 * We do upwind differnecing, so we need to compute the ratio delta_u_upwind /
 * delta_u_local, and how that is computed depends on the local velocity v in
 * the cell. This function is intended for advection purposes, not really for
 * hydro.
 *
 * if v > 0:
 *    compute r = (u_{i} - u_{i-1}) / (u_{i+1} - u_{i})
 * else:
 *    compute r = (u_{i+1} - u_{i+2}) / (u_{i} - u_{i+1})
 *
 * @param pstate* Uip2:  U_{i+2}
 * @param pstate* Uip1:  U_{i+1}
 * @param pstate* Ui:    U_{i}
 * @param pstate* Uim1:  U_{i-1}
 * @param pstate* r:     r for every primitive state
 * @param float vel:     advection velocity
 * ---------------------------------------------------------------------------------------------*/
void limiter_get_r_pstate(pstate *Uip2, pstate *Uip1, pstate *Ui, pstate *Uim1,
                          pstate *r, float vel) {

  if (vel >= 0) {
    r->rho = limiter_r(Ui->rho, Uim1->rho, Uip1->rho);
    r->u[0] = limiter_r(Ui->u[0], Uim1->u[0], Uip1->u[0]);
    r->u[1] = limiter_r(Ui->u[1], Uim1->u[1], Uip1->u[1]);
    r->p = limiter_r(Ui->p, Uim1->p, Uip1->p);
  } else {
    r->rho = limiter_r(Uip1->rho, Uip2->rho, Ui->rho);
    r->u[0] = limiter_r(Uip1->u[0], Uip2->u[0], Ui->u[0]);
    r->u[1] = limiter_r(Uip1->u[1], Uip2->u[1], Ui->u[1]);
    r->p = limiter_r(Uip1->p, Uip2->p, Ui->p);
  }
}


/**
 * Compute the flow parameter r for every component of the conserved states.
 * We always compute r = (u_{i} - u_{i-1}) / (u_{i+1} - u_{i}), and for hydro
 * purposes, we don't really need to care about upwinding.
 *
 * @param cstate* Uip1:  U_{i+1}
 * @param cstate* Ui:    U_{i}
 * @param cstate* Uim1:  U_{i-1}
 * @param cstate* r:     where flow parameter r for every conserved state will
 *                       be stored
 */
void limiter_get_r_cstate(cstate *Uip1, cstate *Ui, cstate *Uim1, cstate *r) {

  r->rho = limiter_r(Ui->rho, Uim1->rho, Uip1->rho);
  r->rhou[0] = limiter_r(Ui->rhou[0], Uim1->rhou[0], Uip1->rhou[0]);
  r->rhou[1] = limiter_r(Ui->rhou[1], Uim1->rhou[1], Uip1->rhou[1]);
  r->E = limiter_r(Ui->E, Uim1->E, Uip1->E);
}


/**
 * in case of advection:
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
 */
float limiter_r(float topleft, float topright, float bottomleft) {

  if (bottomleft == topleft) {
    return ((topleft - topright) * 1e6);
  }
  return ((topleft - topright) / (bottomleft - topleft));
}
