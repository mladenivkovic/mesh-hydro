/* top level file for macro definitions to be used throughout the code*/

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */



#ifndef DEFINES_H
#define DEFINES_H


/* some behaviour options */
/* ---------------------- */

/* keep one global velocity, and keep it constant */
#define ADVECTION_KEEP_VELOCITY_CONSTANT




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
#define GAMMA 1.4

#define GM1 GAMMA-1


/* boundary cells */
#define BC 2 /* how many virtual/ghost boundary cells on each side to make */
#define BCTOT 2*BC


/* boxsize */
#define BOXLEN 1.



/* define solvers as integers */
#define ADVECTION_PWCONST 1
#define ADVECTION_PWLIN 2
#define GODUNOV_UPWIND 3

/* define riemann solvers as integers */
#define NONE 0
#define EXACT 1
#define HLLC 2
#define HLL 3
#define TRRS 4
#define TSRS 5

/* define limiters as integers */
#define MINMOD 1
#define SUPERBEE 2
#define VANLEER 3
#define MC 4

#endif
