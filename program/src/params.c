/* global runtime parameter related stuff. */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#include "params.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

extern params pars;

/**
 * Initialize parameters to default values
 */
void params_init_defaults(void) {

  /* Talking related parameters */
  pars.verbose   = 0;
  pars.nstep_log = 0;

  /* simulation related parameters */
  pars.nsteps = 0;
  pars.tmax   = 0;

  pars.nx       = 100;
  pars.ccfl     = 0.9;
  pars.force_dt = 0;
  pars.boundary = 0;

  pars.nxtot = 100 + BCTOT;
  pars.dx    = BOXLEN / pars.nx;

  /* output related parameters */
  pars.foutput = 0;
  pars.dt_out  = 0;
  strcpy(pars.outputfilename, "");

  strcpy(pars.toutfilename, "");
  pars.use_toutfile = 0;
  pars.noutput_tot  = 0;
  pars.noutput      = 0;
  pars.outputtimes  = NULL;

  /* IC related parameters */
  pars.twostate_ic = 0;
  pars.ndim_ic     = -1;
  strcpy(pars.datafilename, "");

  strcpy(pars.paramfilename, "");

  /* Sources related parameters */
  pars.src_const_acc_x                = 0.;
  pars.src_const_acc_y                = 0.;
  pars.src_const_acc_r                = 0.;
  pars.constant_acceleration          = 0;
  pars.constant_acceleration_computed = 0;
  pars.sources_are_read               = 0;
}


/**
 * Set up parameters that need some processing first
 */
void params_init_derived(void) {

  /* Generate output filename based on ic filename */
  /*-----------------------------------------------*/
  /* do this only if no basename is given */
  if (strlen(pars.outputfilename) == 0) {

    int dot = 0;
    /* extract filename without suffix */
    for (int i = strlen(pars.datafilename); i > 0; i--) {
      if (pars.datafilename[i] == '.') {
        dot = i;
        break;
      }
    }

    int slash = 0;
    /* remove possible directories paths from filename*/
    for (int i = 0; i < (int)strlen(pars.datafilename); i++) {
      if (pars.datafilename[i] == '/') { slash = i; }
    }

    if (dot == 0) dot = strlen(pars.datafilename);
    if (slash > 0) slash += 1;

    char solver[80];
    char riemann[80];
    char limiter[80];

    utils_get_macro_strings(solver, riemann, limiter);

    /* now copy the exact part that you want into filename string */
    strncpy(pars.outputfilename, pars.datafilename + slash, dot - slash);
    pars.outputfilename[dot - slash] = '\0'; /* safety measure */
    strcat(pars.outputfilename, "-");

    /* Add hydro solver */
    strcat(pars.outputfilename, solver);
    strcat(pars.outputfilename, "-");

    /* Add Riemann solver, if present */
    if (strcmp(riemann, "NONE")) { /* 0 if equal */
      strcat(pars.outputfilename, riemann);
      strcat(pars.outputfilename, "-");
    }

    /* Add Limiter */
    if (strcmp(limiter, "NONE")) { /* 0 if equal */
      strcat(pars.outputfilename, limiter);
    } else {
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

  /* -----------------------------------------------------------------------------------
   * If the user set dt_out, pre-compute the output times just like we would
   * have if we used a outputtime file to specify the output times in order to
   * avoid misbehaviour with * floating point divisions
   * -----------------------------------------------------------------------------------
   */

  if (pars.dt_out > 0) {
    /* set the flag */
    pars.use_toutfile = 1;

    /* now get output times */
    int nout;

    if (pars.tmax > 0) {
      nout = floor(pars.tmax / pars.dt_out) + 1; /* add 1 for good measure */
    } else {
      /* assume maximally 9999 outputs */
      nout = 9999;
    }

    pars.outputtimes = malloc(nout * sizeof(float));
    pars.noutput_tot = nout;
    for (int i = 0; i < nout; i++) {
      pars.outputtimes[i] = (i + 1) * pars.dt_out;
    }
  }

  /* Mark if we use constant sources */
  /* ------------------------------- */
#if (SOURCE == SRC_CONST) || (SOURCE == SRC_RADIAL)
  pars.constant_acceleration = 1;
#endif
}


/**
 * Print out current parameters
 */
void params_print_log(void) {

  log_message("\n");
  log_message("Runtime parameters are:\n");
  log_message("\n");

  log_message("Verbose?                     ");
  if (pars.verbose > 0) { printbool(pars.verbose); }
  if (pars.verbose > 0) { printf("\n"); }
  if (pars.nstep_log > 0) {
    log_message("Will write logs every %d steps\n", pars.nstep_log);
  }

  log_message("tmax:                        %g\n", pars.tmax);
  log_message("nsteps:                      %d\n", pars.nsteps);

  log_message("nx:                          %d\n", pars.nx);
  log_message("C_cfl:                       %g\n", pars.ccfl);

  if (pars.force_dt > 0) {
    log_message("Forcing time step size to: %g\n", pars.force_dt);
  }

  log_message("boundary conditions:         ");
  if (pars.verbose > 0) {
    if (pars.boundary == 0) {
      printf("periodic\n");
    } else if (pars.boundary == 1) {
      printf("reflective\n");
    } else if (pars.boundary == 2) {
      printf("transmissive\n");
    }
  }

  log_message("IC file:                     %s\n", pars.datafilename);
  if (pars.twostate_ic) {
    log_message(
      "                             IC file has only two primitive "
      "states.\n"
    );
  }

  if (pars.use_toutfile) {
    if (pars.dt_out == 0.0) {
      log_message("Using output times file:     %s\n", pars.toutfilename);
    } else {
      log_message("dt_out:                      %g\n", pars.dt_out);
    }
  } else {
    log_message("foutput:                     %d\n", pars.foutput);
  }

  log_message("output file basename:        %s\n", pars.outputfilename);

#if SOURCE == SRC_CONST
  log_message("constant source in x:        %g\n", pars.src_const_acc_x);
  log_message("constant source in y:        %g\n", pars.src_const_acc_y);
#elif SOURCE == SRC_RADIAL
  log_message("constant source in r:        %g\n", pars.src_const_acc_r);
#endif

  log_message(
    "----------------------------------------------------------------"
    "-------------------------\n"
  );
}


/**
 * Check whether we can work with these parameters
 */
void params_check(void) {

  log_extra("Checking whether we have valid parameters");

  if (pars.tmax == 0 && pars.nsteps == 0) {
    throw_error(
      "In params_check: I have nsteps = 0 and tmax = 0. You need to "
      "tell me when to stop."
    );
  }

  if (pars.foutput < 0) {
    throw_error("foutput is negative. What do you expect me to do with that?");
  }

  if (pars.dt_out < 0) {
    throw_error("dt_out is negative. What do you expect me to do with that?");
  }

  if (pars.foutput > 0 && pars.dt_out > 0) {
    throw_error(
      "You specified dt_out and foutput > 0. You can't have both, pick one."
    );
  }

  if (pars.use_toutfile) {
    if (pars.dt_out != 0)
      throw_error(
        "You gave me a output time file, but also a dt_out. Decide "
        "which you want and retry."
      );
    if (pars.foutput != 0)
      throw_error(
        "You gave me a output time file, but also a foutput. Decide "
        "which you want and retry."
      );
  }

  if (pars.nx == 0) {
    throw_error("In params_check: I have nx = 0 cells for the sim.");
  }

  if (pars.ndim_ic != NDIM) {
    if (!pars.twostate_ic) {
      int nd = NDIM;
      throw_error(
        "You're trying to use an arbitrary IC filetype for ndim=%d, "
        "but the code is compiled for ndim=%d",
        pars.ndim_ic,
        nd
      );
    }
  }

  /* check source related stuff. */
#ifdef WITH_SOURCES
  if (!pars.sources_are_read) {
    /* Have we compiled with sources, but read in what the sources are? */
    throw_error(
      "Code is compiled to work with sources, but I haven't read in "
      "any source related parameters."
    );
  }
#else
  /* have we compiled without sources, but read in sources? */
  if (pars.sources_are_read) {
    throw_error(
      "Code is compiled to work without sources, but I read in "
      "source related parameters."
    );
  }
#endif
}


/**
 * Check whether we can work with these parameters When using the program as a
 * riemann solver
 */
void params_check_riemann(void) {

  log_extra("Checking whether we have valid parameters");

  if (pars.tmax == 0) {
    throw_error(
      "You need to specify tmax so I know at what time to sample the "
      "solution."
    );
  }

  if (pars.nx == 0) {
    throw_error(
      "In params_check: I have nx = 0 cells for the sim. You need to "
      "tell me how many cells you want."
    );
  }
}


/**
 * Create a filename for the output for when the code is employed as a Riemann
 * solver only
 */
void params_generate_riemann_output_filename(void) {

  /* do this only if no basename is given */
  if (strlen(pars.outputfilename) == 0) {

    int dot = 0;
    /* extract filename without suffix */
    for (int i = strlen(pars.datafilename); i > 0; i--) {
      if (pars.datafilename[i] == '.') {
        dot = i;
        break;
      }
    }

    int slash = 0;
    /* remove possible directories paths from filename*/
    for (int i = 0; i < (int)strlen(pars.datafilename); i++) {
      if (pars.datafilename[i] == '/') { slash = i; }
    }

    if (dot == 0) { dot = strlen(pars.datafilename); }
    if (slash > 0) { slash += 1; }

    char solver[80];
    char riemann[80];
    char limiter[80];

    utils_get_macro_strings(solver, riemann, limiter);

    /* now copy the exact part that you want into filename string */
    strcpy(pars.outputfilename, "");
    strncpy(pars.outputfilename, pars.datafilename + slash, dot - slash);
    pars.outputfilename[dot - slash] = '\0';
    strcat(pars.outputfilename, "-RIEMANN-");
    strcat(pars.outputfilename, riemann);
  }
}
