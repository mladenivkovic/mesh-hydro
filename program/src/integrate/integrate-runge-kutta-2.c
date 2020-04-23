/* Runge-Kutta 2 integrator             */

/* Written by Mladen Ivkovic, APR 2020
 * mladen.ivkovic@hotmail.com           */


#include "cell.h"
#include "gas.h"
#include "integrate.h"
#include "params.h"
#include "sources.h"



void integrate(cstate* U, float acc[2], float dt, cstate* Unew){
  /* ------------------------------------------
   * integrate dU/dt = S over the time step dt
   * ------------------------------------------ */

  /* get source vector */
  cstate S;
  gas_init_cstate(&S);
  sources_get_source_vector(&S, acc, U);

  /* compute first midstep */
  cstate k1;
  k1.rho     = S.rho     * dt;
  k1.rhou[0] = S.rhou[0] * dt;
  k1.rhou[1] = S.rhou[1] * dt;
  k1.E       = S.E       * dt;

  /* get intermediate state to be able to compute next step */
  cstate Umid;
  Umid.rho     = U->rho     + k1.rho;
  Umid.rhou[0] = U->rhou[0] + k1.rhou[0];
  Umid.rhou[1] = U->rhou[1] + k1.rhou[1];
  Umid.E       = U->E       + k1.E;

  /* get intermediate step Source vector */
  sources_get_source_vector(&S, acc, &Umid);

  /* get second intermediate step */
  cstate k2;
  k2.rho     = S.rho     * dt;
  k2.rhou[0] = S.rhou[0] * dt;
  k2.rhou[1] = S.rhou[1] * dt;
  k2.E       = S.E       * dt;


  /* get solution */
  Unew->rho     = U->rho     + 0.5 * (k1.rho     + k2.rho);
  Unew->rhou[0] = U->rhou[0] + 0.5 * (k1.rhou[0] + k2.rhou[0]);
  Unew->rhou[1] = U->rhou[1] + 0.5 * (k1.rhou[1] + k2.rhou[1]);
  Unew->E       = U->E       + 0.5 * (k1.E       + k2.E);

}
