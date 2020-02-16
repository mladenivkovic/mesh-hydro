/* Relations, states etc concerning ideal gasses */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */


#ifndef GAS_H
#define GAS_H

/* primitive state */
typedef struct {
  float rho;  /* density */
  float ux;   /* velocity in x direction */
  float uy;   /* velocity in y direction */
  float p;    /* pressure */
} pstate;


/* conserved state */
typedef struct {
  float rho;    /* density */
  float rhoux;  /* specific momentum in x direction */
  float rhouy;  /* specific momentum in x direction */
  float E;      /* specific energy */
} cstate; 



void gas_init_pstate(pstate *p);
void gas_init_cstate(cstate *c);

float gas_soundspeed(pstate* s);
float gas_energy(pstate* s);


#endif
