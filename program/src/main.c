/* =======================================================
 * mesh based Conservation Law solver !
 * run with ./hydro parameterfile.txt ICfile.dat
 * See the README.md in parent directory on how the
 * parameter and IC file can look like.
 *
 * Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com
 * ======================================================= */

#include <stdio.h>
#include <time.h> /* measure time */

#include "cell.h"
#include "io.h"
#include "params.h"
#include "solver.h"
#include "utils.h"

/* ------------------ */
/* Initialize globals */
/* ------------------ */

params pars;

#if NDIM == 1
cell* grid;
#elif NDIM == 2
cell** grid;
#endif


int main(int argc, char* argv[]) {

  /* timing stuff */
  clock_t step_start, step_end;
  clock_t all_start, all_end;

  all_start = clock();

  /* Unnecessary things first :) */
  print_header();

  /* initialize parameters struct */
  params_init_defaults();

  /* read cmdline args and parameter file */
  io_read_cmdlineargs(argc, argv);
  io_read_paramfile();

  /* how many lines to skip next time you open the IC file */
  int skiplines_ic = 0;
  /* check which IC type we have */
  io_read_ic_type(&skiplines_ic);

  /* read in output times if necessary */
  if (pars.use_toutfile) { io_read_toutfile(); }

  params_check();        /* check whether we can work with this setup. */
  params_init_derived(); /* process the parameters you got. */

  /* print / announce stuff for logging */
  print_compile_defines();
  params_print_log();

  /* initialize the grid */
  cell_init_grid();

  /* read in the full IC file */
  if (pars.twostate_ic) {
    io_read_ic_twostate(skiplines_ic);
  } else {
    io_read_ic_arbitrary(skiplines_ic);
  }

#ifdef ADVECTION_KEEP_VELOCITY_CONSTANT
  solver_advection_check_global_velocity();
#endif
  /* translate the read-in primitive vars to conservative ones */
  cell_get_cstates_from_pstates();

  /* Initialize counters and time */
  int   step     = 0; /* step counter */
  int   outcount = 0; /* number of the output that we're writing */
  float t        = 0; /* time */
  float dt       = 0; /* time step size */

  int write_output = 0; /* whether the time step was reduced because we need to
                           write an output */

  float mtot_init = cell_get_total_mass(); /* for checks every step */

  log_extra("Writing initial output");
  io_write_output(&outcount, step, t);

  log_message("\n");
  log_message(
    "%14s %14s %14s %14s  %14s\n",
    "step",
    "time",
    "dt",
    "m_now/m_ini",
    "time step took"
  );

  /* --------------------
   *   Main loop
   * -------------------- */
  while (1) {
    if (pars.tmax > 0 && t >= pars.tmax) break;
    if (pars.nsteps > 0 && step == pars.nsteps) break;

    step_start = clock(); /* timer */

    /* where the actual magic happens */
    solver_step(&t, &dt, step, &write_output);

    step_end = clock(); /* timer */

    /* and update time and step*/
    t += dt;
    step += 1;

    /* write output if you have to */
    if (write_output) { io_write_output(&outcount, step, t); }

    /* announce */
    if (pars.nstep_log == 0 || step % pars.nstep_log == 0) {
      log_message(
        "%14d %14.6e %14.6e %14.6e %14.3es\n",
        step,
        t,
        dt,
        cell_get_total_mass() / mtot_init,
        (float)(step_end - step_start) / CLOCKS_PER_SEC
      );
    }
  }

  /* if you haven't written the output in the final step, do it now */
  if (!write_output) { io_write_output(&outcount, step, t); }

  all_end = clock();

  /* don't use log_message, I want final stats even for verbose = 0 */
  printf("\n");
  printf("  Finished clean. Yay!\n");
  printf("  Final stats:\n");
  printf("\n");
  printf(
    "    Total runtime was       %12.6fs\n",
    (float)(all_end - all_start) / CLOCKS_PER_SEC
  );
  printf(
    "    m_now/m_ini =           %12.6f\n", cell_get_total_mass() / mtot_init
  );
  printf("    final number of steps = %12d\n", step);

  return (0);
}
