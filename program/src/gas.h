/* Relations, states etc concerning ideal gasses */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */


#ifndef GAS_H
#define GAS_H

/* primitive state */
typedef struct {
  float rho;  /* density */
  float* u;   /* velocity vector. u[0] = ux, u[1] = uy */
  float p;    /* pressure */
} pstate;


/* conserved state */
typedef struct {
  float rho;    /* density */
  float *rhou;  /* specific momentum. rhou[0] in x direction, rhou[1] in y direction */
  float E;      /* specific energy */
} cstate; 



void gas_init_pstate(pstate *p);
void gas_init_cstate(cstate *c);

float gas_soundspeed(pstate *s);
float gas_energy(pstate* s);


#endif
