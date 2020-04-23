/* top level file for macro definitions to be used throughout the code*/

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */



#ifndef DEFINES_H
#define DEFINES_H


/* -------------------------------------------------------------------
 * some behaviour options you are free to tinker with
 * ------------------------------------------------------------------- */

/* keep one global velocity, and keep it constant */
#ifdef ADVECTION
#define ADVECTION_KEEP_VELOCITY_CONSTANT
#endif

/* whether to use an adaptive, or a simple wave speed estimate for the HLLC
 * Riemann solver */
#define HLLC_USE_ADAPTIVE_SPEED_ESTIMATE


/* Physical constants */

#define GAMMA  (5./3.)
// #define GAMMA 1.4

/* boundary cells */
#define BC 2 /* how many virtual/ghost boundary cells on each side to make */
#define BCTOT 2*BC

/* boxsize */
#define BOXLEN 1.

/* define slope of each cell in MUSCL scheme as
 * slope = 0.5* (1 + omega) (U_{i}-U_{i-1}) +  0.5 * (1 - omega) (U_{i+1} - U_{i}) */
#define OMEGA 0.




/* ---------------------------------------------------------------------------------- *
 * you shouldn't be modifying stuff below
 * this line unless you know what you're doing
 * ---------------------------------------------------------------------------------- */

/* File related stuff */

#define MAX_FNAME_SIZE 200    /* limit for file name size */
#define MAX_LINE_SIZE 200     /* limit for line length in formatted file which is read in*/


/* Macro functions */

#define STR(x) STR_(x)
#define STR_(x) #x

/* iteration tolerance */
#define EPSILON_ITER 1e-6


/* minimal timestep size */
#define DT_MIN 1e-10

/* "nonzero zeros" for vacuum */
#ifdef USE_AS_RIEMANN_SOLVER
/* set the "small" values to actually zero, so that only correct vacuum sates 
 * are recognized as such */
#define SMALLRHO 0.
#define SMALLU 0.
#define SMALLP 0.
#else
/* cheat for stability in Godunov type finite volume schemes*/
#define SMALLRHO 1e-6
#define SMALLU 1e-6
#define SMALLP 1e-6
#endif








/* ---------------------------------------------------------------------------------- *
 * Nobody should be changing things below this line
 * ----------------------------------------------------------------------------------*/


/* Compute gamma related constants */

static const float GM1 = GAMMA - 1.;
static const float GP1 = GAMMA + 1.;
static const float GP1OGM1 = (GAMMA + 1.) / (GAMMA - 1.);
static const float GM1OGP1 = (GAMMA - 1.) / (GAMMA + 1.);
static const float ONEOVERGAMMA = 1. / GAMMA;
static const float GM1HALF = 0.5 * (GAMMA - 1.);
static const float BETA = 0.5 * (GAMMA - 1.) / GAMMA;


/* PI */
#define PI 3.14159265 /* I know there are better ways of doing this. Sue me. */

/* define solvers as integers */
#define ADVECTION_PWCONST 1
#define ADVECTION_PWLIN 2
#define GODUNOV 3
#define ADVECTION_WAF 4
#define WAF 5
#define MUSCL 6

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

/* define sources as integers */
#define SRC_CONST 1
#define SRC_RADIAL 2


#endif
