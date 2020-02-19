/* Relations, states etc concerning ideal gasses */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */


#ifndef GAS_H
#define GAS_H

#include "defines.h"

/* primitive state */
typedef struct {
  MYFLOAT rho;  /* density */
  MYFLOAT ux;   /* velocity in x direction */
  MYFLOAT uy;   /* velocity in y direction */
  MYFLOAT p;    /* pressure */
} pstate;


/* conserved state */
typedef struct {
  MYFLOAT rho;    /* density */
  MYFLOAT rhoux;  /* specific momentum in x direction */
  MYFLOAT rhouy;  /* specific momentum in x direction */
  MYFLOAT E;      /* specific energy */
} cstate; 



void gas_init_pstate(pstate *p);
void gas_init_cstate(cstate *c);

MYFLOAT gas_soundspeed(pstate* s);
MYFLOAT gas_energy(pstate* s);


#endif
