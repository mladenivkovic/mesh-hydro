/* top level file for macro definitions to be used throughout the code*/

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */



#ifndef DEFINES_H
#define DEFINES_H


/* some behaviour options */
/* ---------------------- */

/* keep one global velocity, and keep it constant */
#ifdef ADVECTION
#define ADVECTION_KEEP_VELOCITY_CONSTANT
#endif




/* -------------------------------------- */
/* you shouldn't be modifying stuff below
 * this line unless you know what you're
 * doing                                  */

/* File related stuff */

#define MAX_FNAME_SIZE 200    /* limit for file name size */
#define MAX_LINE_SIZE 200     /* limit for line length in formatted file which is read in*/


/* Macro functions */

#define STR(x) STR_(x)
#define STR_(x) #x


/* Physical constants */

#define PI 3.14159265358979
#define GAMMA (5./3.)

#define GM1 (GAMMA-1.)
#define GP1 (GAMMA+1.)
#define GP1OGM1 ((GP1)/(GM1))
#define GM1OGP1 ((GM1)/(GP1))
#define ONEOVERGAMMA (1./GAMMA)
#define GM1HALF (0.5*(GM1))
#define BETA ((GM1HALF) / GAMMA)


/* boundary cells */
#define BC 2 /* how many virtual/ghost boundary cells on each side to make */
#define BCTOT 2*BC


/* boxsize */
#define BOXLEN 1.


/* iteration tolerance */
#define EPSILON_ITER 1e-6


/* define solvers as integers */
#define ADVECTION_PWCONST 1
#define ADVECTION_PWLIN 2
#define GODUNOV 3

/* define riemann solvers as integers */
#define NONE 0
#define EXACT 1
#define HLLC 2
#define TRRS 3
#define TSRS 4

/* define limiters as integers */
#define MINMOD 1
#define SUPERBEE 2
#define VANLEER 3
#define MC 4

#endif
