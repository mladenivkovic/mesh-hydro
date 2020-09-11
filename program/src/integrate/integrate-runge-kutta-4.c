/* Runge-Kutta 4 integrator             */

/* Written by Mladen Ivkovic, APR 2020
 * mladen.ivkovic@hotmail.com           */

#include "cell.h"
#include "gas.h"
#include "integrate.h"
#include "params.h"
#include "sources.h"

void integrate(cstate *U, float acc[2], float dt, cstate *Unew) {
  /* ------------------------------------------
   * integrate dU/dt = S over the time step dt
   * ------------------------------------------ */

  /* get source vector */
  cstate S;
  gas_init_cstate(&S);
  sources_get_source_vector(&S, acc, U);

  /* compute first midstep */
  cstate k1;
  k1.rho = S.rho * dt;
  k1.rhou[0] = S.rhou[0] * dt;
  k1.rhou[1] = S.rhou[1] * dt;
  k1.E = S.E * dt;

  /* get intermediate state to be able to compute next step */
  cstate Umid;
  Umid.rho = U->rho + 0.5 * k1.rho;
  Umid.rhou[0] = U->rhou[0] + 0.5 * k1.rhou[0];
  Umid.rhou[1] = U->rhou[1] + 0.5 * k1.rhou[1];
  Umid.E = U->E + 0.5 * k1.E;

  /* get intermediate step Source vector */
  sources_get_source_vector(&S, acc, &Umid); /* should be at t+0.5dt */

  /* get second intermediate step k2 */
  cstate k2;
  k2.rho = S.rho * dt;
  k2.rhou[0] = S.rhou[0] * dt;
  k2.rhou[1] = S.rhou[1] * dt;
  k2.E = S.E * dt;

  /* get intermediate state to be able to compute next step */
  Umid.rho = U->rho + 0.5 * k2.rho;
  Umid.rhou[0] = U->rhou[0] + 0.5 * k2.rhou[0];
  Umid.rhou[1] = U->rhou[1] + 0.5 * k2.rhou[1];
  Umid.E = U->E + 0.5 * k2.E;

  /* get intermediate step Source vector */
  sources_get_source_vector(&S, acc, &Umid); /* should be at t+0.5dt */

  /* get second intermediate step k2 */
  cstate k3;
  k3.rho = S.rho * dt;
  k3.rhou[0] = S.rhou[0] * dt;
  k3.rhou[1] = S.rhou[1] * dt;
  k3.E = S.E * dt;

  /* get intermediate state to be able to compute next step */
  Umid.rho = U->rho + k3.rho;
  Umid.rhou[0] = U->rhou[0] + k3.rhou[0];
  Umid.rhou[1] = U->rhou[1] + k3.rhou[1];
  Umid.E = U->E + k3.E;

  /* get intermediate step Source vector */
  sources_get_source_vector(&S, acc, &Umid); /* should be at t+dt */

  /* get second intermediate step k2 */
  cstate k4;
  k4.rho = S.rho * dt;
  k4.rhou[0] = S.rhou[0] * dt;
  k4.rhou[1] = S.rhou[1] * dt;
  k4.E = S.E * dt;

  /* get solution */
  Unew->rho = U->rho + (k1.rho + 2 * k2.rho + 2 * k3.rho + k4.rho) / 6;
  Unew->rhou[0] =
      U->rhou[0] +
      (k1.rhou[0] + 2 * k2.rhou[0] + 2 * k3.rhou[0] + k4.rhou[0]) / 6;
  Unew->rhou[1] =
      U->rhou[1] +
      (k1.rhou[1] + 2 * k2.rhou[1] + 2 * k3.rhou[1] + k4.rhou[1]) / 6;
  Unew->E = U->E + (k1.E + 2 * k2.E + 2 * k3.E + k4.E) / 6;
}
