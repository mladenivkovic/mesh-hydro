/* no limiter  */


/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */

#include "gas.h"
#include "cell.h"
#include "limiter.h"
#include "params.h"


#if NDIM == 1
extern cell *grid;
#elif NDIM == 2
extern cell **grid;
#endif


extern params pars;


void limiter_get_slope_left(cell* c, pstate* slope, int dimension){
  /* ------------------------------------------------------------------------
   * Compute the left slope of given cell c, i.e. the slope for the flux
   * F_{i-1/2}.
   * ------------------------------------------------------------------------ */

  int i, j;
  cell* target;
  cell_get_ij(c, &i, &j);

  if (dimension == 0){
#if NDIM == 1
    if (c->prim.u[0] >= 0){
      target = &grid[i-1];
    } else {
      target = &grid[i];
    }
#elif NDIM == 2
    if (c->prim.u[0] >= 0){
      target = &grid[i-1][j];
    } else {
      target = &grid[i][j];
    }
  } else if (dimension == 1){
    if (c->prim.u[1] >= 0){
      target = &grid[i][j-1];
    } else {
      target = &grid[i][j];
    }
#endif
  }

  limiter_get_slope(target, slope, dimension);
}



void limiter_get_slope_right(cell* c, pstate* slope, int dimension){
  /* ------------------------------------------------------------------------
   * Compute the right slope of given cell c, i.e. the slope for the flux
   * F_{i+1/2}.
   * ------------------------------------------------------------------------ */

  int i, j;
  cell* target;
  cell_get_ij(c, &i, &j);

  if (dimension == 0){
#if NDIM == 1
    if (c->prim.u[0] >= 0){
      target = &grid[i];
    } else {
      target = &grid[i+1];
    }
#elif NDIM == 2
    if (c->prim.u[0] >= 0){
      target = &grid[i][j];
    } else {
      target = &grid[i+1][j];
    }
  } else if (dimension == 1){
    if (c->prim.u[1] >= 0){
      target = &grid[i][j];
    } else {
      target = &grid[i][j+1];
    }
#endif
  }

  limiter_get_slope(target, slope, dimension);
}




void limiter_get_slope(cell* c, pstate* slope, int dimension){
  /* -----------------------------------------------------------
   * compute the slope of cell c and write it in the slope.
   * dimension: 0 = x, 1 = y
   * Here, we just use a centered difference
   * ----------------------------------------------------------- */

  int i, j;
  cell left, right;

  cell_get_ij(c, &i, &j);

  /* first find the cells left and right of you */
  if (dimension == 0){

#if NDIM == 1
    right = grid[i+1];
    left = grid[i-1];
#elif NDIM == 2
    right = grid[i+1][j];
    left = grid[i-1][j];
#endif

  } else if (dimension == 1){

#if NDIM == 2
    right = grid[i][j+1];
    left = grid[i][j-1];
#endif
  }

  /* and now do what shall be done */
  slope->rho = 0.5 * (right.prim.rho - left.prim.rho) / pars.dx;
  slope->u[0] = 0.5 * (right.prim.u[0] - left.prim.u[0]) / pars.dx;
  slope->u[1] = 0.5 * (right.prim.u[1] - left.prim.u[1]) / pars.dx;
  slope->p = 0.5 * (right.prim.p - left.prim.p) / pars.dx;
}
