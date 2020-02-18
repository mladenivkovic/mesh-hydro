/* All around parameters used in the simulation. */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */


#ifndef PARAMS_H
#define PARAMS_H

#include "defines.h"

/* GLOBAL PARAMETERS */
typedef struct {
  int verbose;                          /* how talkative I am */

  int nsteps;                           /* how many steps to take */
  float tmax;                           /* at what time to end */

  float ccfl;                           /* CFL coefficient */
  int nx;                               /* number of mesh points */
  int nxtot;                            /* number of mesh points, including boundary cells */
  float dx;                             /* cell size */
  float force_dt;                       /* force a time step size (except if you need to write an output) */
  
  int twostate_ic;                      /* whether IC are left/right state only */
  int ndim_ic;                          /* dimension of IC file */
  char datafilename[MAX_FNAME_SIZE];    /* IC data filename */

  int foutput;                          /* after how many steps to write output */
  float dt_out;                         /* time interval between outputs */
  char outputfilename[MAX_FNAME_SIZE];  /* Output file name basename */
  int nstep_log;                        /* on how many steps to write log message of time step */

  char toutfilename[MAX_FNAME_SIZE];    /* file name containing output times */
  int use_toutfile;                     /* whether we're using the t_out_file */
  int noutput_tot;                      /* how many outputs we will be writing. Only used if(use_toutfile) */
  int noutput;                          /* at which output we are. Only used if(use_toutfile) */
  float *outputtimes;                   /* array of output times given in the output file */

  char paramfilename[MAX_FNAME_SIZE];   /* parameter filename */

  int boundary;                         /* boundary condition */
} params;






void params_init_defaults();
void params_init_derived();
void params_print_log();
void params_check();


#endif
