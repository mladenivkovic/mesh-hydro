/* global runtime parameter related stuff. */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "params.h"
#include "utils.h"


extern params pars;





void params_init_defaults(){
  /*------------------------------------------*/
  /* Initialize parameters to default values  */
  /*------------------------------------------*/

  pars.verbose = 0;

  pars.nsteps = 1;
  pars.tmax = 0;

  pars.ccfl = 0.9;
  pars.nx = 100;
  pars.nxtot = 100+BCTOT;
  pars.dx = BOXLEN / 100;

  pars.twostate_ic = 0;
  pars.ndim_ic = -1;
  strcpy(pars.datafilename, "");

  pars.foutput = 0;
  pars.dt_out = 0;
  strcpy(pars.outputfilename, "");
  pars.nstep_log = 0;

  pars.use_toutfile = 0;
  pars.noutput_tot = 0;
  pars.noutput = 0;
  strcpy(pars.toutfilename, "");
  pars.outputtimes = NULL;

  strcpy(pars.paramfilename, "");

  pars.boundary = 0; 
}



void params_init_derived(){
  /* ------------------------------------------------- */
  /* Set up parameters that need some processing first */
  /* ------------------------------------------------- */

  /* Generate output filename based on ic filename */
  /*-----------------------------------------------*/
  if (strlen(pars.outputfilename)==0) {

    int dot = 0;
    /* extract filename without suffix */
    for (int i = strlen(pars.datafilename); i > 0; i--){
      if (pars.datafilename[i] == '.'){
        dot = i;
        break;
      }
    }

    int slash = 0;
    /* remove possible directories paths from filename*/
    for (int i = 0; i < (int) strlen(pars.datafilename); i++){
      if (pars.datafilename[i] == '/'){
        slash = i;
      }
    }

    if (dot==0) dot = strlen(pars.datafilename);
    if (slash > 0) slash += 1;

    char solver[80];
    char riemann[80];
    char limiter[80];

    utils_get_macro_strings(solver, riemann, limiter);

    /* now copy the exact part that you want into filename string */
    strncpy(pars.outputfilename, pars.datafilename+slash, dot-slash);
    strcat(pars.outputfilename, "-");

    /* Add hydro solver */
    strcat(pars.outputfilename, solver);
    strcat(pars.outputfilename, "-");

    /* Add Riemann solver, if present */
    if(strcmp(riemann, "NONE")){ /* 0 if equal */
      strcat(pars.outputfilename, riemann);
      strcat(pars.outputfilename, "-");
    }

    /* Add Limiter */
    if(strcmp(limiter, "NONE")){ /* 0 if equal */
      strcat(pars.outputfilename, limiter);
    }
    else {
      strcat(pars.outputfilename, "NO_LIMITER");
    }
    strcat(pars.outputfilename, "-");

    /* Add dimension */
    strcat(pars.outputfilename, STR(NDIM));
    strcat(pars.outputfilename, "D");



  }

  /* Compute total number of cells per dimension */
  pars.nxtot = pars.nx + BCTOT;

  /* Compute dx */
  pars.dx = BOXLEN / pars.nx;
}




void params_print_log(){
  /*------------------------------------------*/
  /* Print out current parameters             */
  /*------------------------------------------*/

  log_message("\n");
  log_message("Runtime parameters are:\n");
  log_message("\n");

  log_message("Verbose?                     ");
  if (pars.verbose > 0) printbool(pars.verbose);
  if (pars.verbose > 0) printf("\n");

  log_message("tmax:                        %g\n", pars.tmax);
  log_message("nsteps:                      %d\n", pars.nsteps);

  log_message("nx:                          %d\n", pars.nx);
  log_message("C_cfl:                       %g\n", pars.ccfl);

  log_message("boundary conditions:         ");
  if (pars.verbose > 0) {
    if (pars.boundary == 0){
      printf("periodic\n");
    } else if (pars.boundary == 1){
      printf("reflective\n");
    } else if (pars.boundary == 2){
      printf("transmissive\n");
    }
  }

  log_message("IC file:                     %s\n", pars.datafilename);
  if (pars.twostate_ic){
  log_message("                             IC file has only two primitive states.\n");
  }

  if (pars.use_toutfile){
    log_message("Using output times file:     %s\n", pars.toutfilename);
  } else{
    if (pars.dt_out == 0.0){
      log_message("foutput:                     %d\n", pars.foutput);
    } else {
      log_message("dt_out:                      %d\n", pars.dt_out);
    }
  }
  if (pars.nstep_log > 0){
    log_message("Will write logs every %d steps\n", pars.nstep_log);
  }
  log_message("output file basename:        %s\n", pars.outputfilename);
  log_message("-----------------------------------------------------------------------------------------\n");
}




void params_check(){
  /* ----------------------------------------------- */
  /* Check whether we can work with these parameters */
  /* ----------------------------------------------- */

  log_extra("Checking whether we have valid parameters");

  if (pars.tmax == 0 && pars.nsteps == 0){
    throw_error("In params_check: I have nsteps = 0 and tmax = 0. You need to tell me when to stop.");
  }

  if (pars.foutput < 0) {
    throw_error("foutput is negative. What do you expect me to do with that?");
  }

  if (pars.dt_out < 0) {
    throw_error("dt_out is negative. What do you expect me to do with that?");
  }

  if (pars.foutput > 0 && pars.dt_out > 0) {
    throw_error("You specified dt_out and foutput > 0. You can't have both, pick one.");
  }

  if (pars.use_toutfile){
    if (pars.dt_out != 0) throw_error("You gave me a output time file, but also a dt_out. Decide which you want and retry.");
    if (pars.foutput != 0) throw_error("You gave me a output time file, but also a foutput. Decide which you want and retry.");
  }


  if (pars.nx == 0) {
    throw_error("In params_check: I have nx = 0 cells for the sim.");
  }

  if (pars.ndim_ic != NDIM){
    if (!pars.twostate_ic) {
      int nd = NDIM;
      throw_error("You're trying to use an arbitrary IC filetype for ndim=%d, but the code is compiled for ndim=%d", pars.ndim_ic, nd);
    }
  }
}
