/* Relations, states etc concerning ideal gasses */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef GAS_H
#define GAS_H

/* primitive state */
typedef struct {
  float rho;  /* density */
  float u[2]; /* velocity vector. u[0] = ux, u[1] = uy */
  float p;    /* pressure */
} pstate;

/* conserved state */
typedef struct {
  float rho;     /* density */
  float rhou[2]; /* specific momentum. rhou[0] in x direction, rhou[1] in y
                    direction */
  float E;       /* specific energy */
} cstate;

void gas_init_pstate(pstate *p);
void gas_init_cstate(cstate *c);
void gas_prim_to_cons(pstate *p, cstate *c);
void gas_cons_to_prim(cstate *c, pstate *p);
void gas_get_cflux_from_pstate(pstate *p, cstate *f, int dimension);
void gas_get_cflux_from_cstate(cstate *c, cstate *f, int dimension);

float gas_soundspeed(pstate *s);
float gas_energy(pstate *s);

#endif
