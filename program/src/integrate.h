/* top level file for integrators */

/* Written by Mladen Ivkovic, APR 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

void integrate(cstate *U, float acc[2], float dt, cstate *Unew);

#endif
