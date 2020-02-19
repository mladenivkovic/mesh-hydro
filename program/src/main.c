/* =======================================================
 * mesh based Conservation Law solver !
 * run with ./hydro parameterfile.txt ICfile.dat
 * See the README.md in parent directory on how the
 * parameter and IC file can look like.
 *
 * Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com
 * ======================================================= */


#include <stdlib.h>
#include <stdio.h>
#include <time.h>  /* measure time */

#include "cell.h"
#include "defines.h"
#include "gas.h"
#include "solver.h"
#include "io.h"
#include "params.h"
#include "riemann.h"
#include "utils.h"




/* ------------------ */
/* Initialize globals */
/* ------------------ */

params pars;

#if NDIM == 1
cell *grid;
#elif NDIM == 2
cell **grid;
#endif






/* ====================================== */
int main(int argc, char* argv[]){
/* ====================================== */

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

  /* check which IC type we have */
  int skiplines_ic = 0; /* how many lines to skip next time you open the IC file */
  io_read_ic_type(&skiplines_ic);

  /* read in output times if necessary */
  if (pars.use_toutfile) io_read_toutfile();

  params_init_derived(); /* process the parameters you got. */
  params_check();        /* check whether we can work with this setup. */

  /* print / announce stuff for logging */
  print_compile_defines();
  params_print_log();



  /* initialize the grid */
  cell_init_grid();

  /* read in the full IC file */
  if (pars.twostate_ic){
    io_read_ic_twostate(skiplines_ic);
  } else {
    io_read_ic_arbitrary(skiplines_ic);
  }

#ifdef ADVECTION_KEEP_VELOCITY_CONSTANT
  solver_advection_check_global_velocity();
#endif



  /* Initialize counters and time */
  int step = 0;         /* step counter */
  int outcount = 0;     /* number of the output that we're writing */
  MYFLOAT t = 0;          /* time */
  MYFLOAT dt = 0;         /* time step size */

  int write_output = 0; /* whether the time step was reduced because we need to write an output */

  MYFLOAT mtot_init = cell_get_total_mass(); /* for checks every step */

  log_extra("Writing initial output");
  io_write_output(&outcount, step, t);


  /* -------------------- 
   *   Main loop
   * -------------------- */
  while(1) {
    if (pars.tmax > 0 && t >= pars.tmax) break; 
    if (pars.nsteps>0 && step == pars.nsteps) break;


    step_start = clock(); /* timer */

    /* where the actual magic happens */
    solver_step(&t, &dt, step,  &write_output);

    step_end = clock(); /* timer */

    /* and update time and step*/
    t += dt;
    step += 1;

    /* write output if you have to */
    if (write_output){
      io_write_output(&outcount, step, t);
    }

    /* announce */
    if (pars.nstep_log == 0 || step % pars.nstep_log == 0) {
      log_message("%10d t = %12.6f dt = %12.6f m_tot = %12.6f step took %12.6fs\n",
                  step, t, dt, cell_get_total_mass()/mtot_init, 
                  (MYFLOAT)(step_end - step_start) / CLOCKS_PER_SEC);
    }
  }

  /* if you haven't written the output in the final step, do it now */
  if (!write_output){ io_write_output(&outcount, step, t); }


  all_end = clock();

  printf("\n");
  printf("  Finished clean. Yay!\n");
  printf("  Total runtime was %12.6fs, mtot = %12.6f, nsteps = %12d\n", 
              (MYFLOAT)(all_end - all_start)/CLOCKS_PER_SEC, cell_get_total_mass()/mtot_init, step);

  return(0);
}
