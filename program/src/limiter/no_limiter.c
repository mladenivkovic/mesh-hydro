/* no limiter  */


/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */

#include "gas.h"
#include "cell.h"
#include "limiter.h"
#include "params.h"
#include "utils.h"


#if NDIM == 1
extern cell *grid;
#elif NDIM == 2
extern cell **grid;
#endif


extern params pars;




void limiter_get_psi_no_limiter(cell* c, pstate* psi, float cfl, int dimension){
  /* ------------------------------------------------------------------------
   * Compute the flux limiter psi_{i+1/2}
   *
   * cell* c:     for which cell i to work for
   * pstate* psi: where the flux limiter will be stored
   * clf:         a * dt / dx
   * dimension:   for which dimension we're working
   * ------------------------------------------------------------------------ */

  psi->rho = 1.;
  psi->u[0] = 1.;
  psi->u[1] = 1.;
  psi->p = 1.;
}





float limiter_psi_of_r(float r){
  /*-------------------------------------------
   * compute the actual psi(r) 
   *
   * NOTE: This should actually
   * never be called for LIMITER == NONE
   * ------------------------------------------*/
  throw_error("You shouldn't be calling limiter_psi_of_r for LIMITER = NONE");
  return(1.);
}
