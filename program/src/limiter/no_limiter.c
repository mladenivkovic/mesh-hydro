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




void limiter_get_phi_no_limiter(cell* c, pstate* phi, int dimension){
  /* ------------------------------------------------------------------------
   * Compute the flux limiter phi_{i+1/2}
   *
   * cell* c:     for which cell i to work for
   * pstate* phi: where the flux limiter will be stored
   * dimension:   for which dimension we're working
   * ------------------------------------------------------------------------ */

  phi->rho = 1.;
  phi->u[0] = 1.;
  phi->u[1] = 1.;
  phi->p = 1.;
}





float limiter_phi_of_r(float r){
  /*-------------------------------------------
   * compute the actual phi(r) 
   *
   * NOTE: This should actually
   * never be called for LIMITER == NONE
   * ------------------------------------------*/
  throw_error("You shouldn't be calling limiter_phi_of_r for LIMITER = NONE");
  return(1.);
}
