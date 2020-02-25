/* =======================================================
 * Use the code to only solve a Riemann problem.
 * run with ./riemann parameterfile.txt ICfile.dat
 * See the README.md in parent directory on how the
 * parameter and IC file can look like.
 *
 * Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com
 * ======================================================= */


#include <stdlib.h>
#include <stdio.h>

#include "cell.h"
#include "defines.h"
#include "gas.h"
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

  /* initialize parameters struct */
  params_init_defaults();

  /* read cmdline args and parameter file */
  io_read_cmdlineargs(argc, argv);
  io_read_paramfile();

  /* check which IC type we have */
  int skiplines_ic = 0; /* how many lines to skip next time you open the IC file */
  io_read_ic_type(&skiplines_ic);

  if (!pars.twostate_ic){
    printf("I can only solve twostate-ICs.");
    exit(2);
  }

  params_init_derived(); /* process the parameters you got. */
  params_check_riemann(); /* check whether we can work with this setup. */


  /* print / announce stuff for logging */
  print_compile_defines();
  params_print_log();



  /* initialize the grid */
  cell_init_grid();

  /* read in the full IC file */
  io_read_ic_twostate(skiplines_ic);

  /* get the left and right state of the original problem */
  pstate left = grid[pars.nx/2+BC-1].prim;
  pstate right = grid[pars.nx/2+BC].prim;

  /* pretend we're doing hydro */
  int outcount = 0;
  int step = 0;
  float t = 0;

  log_extra("Writing initial output");
  io_write_output(&outcount, step, t);



  outcount = 1;
  step = 1;
  t = pars.tmax;


  log_message("Solving Riemann problem.\n");
  log_message("rho_L %12.6f; rho_R %12.6f\n", left.rho, right.rho);
  log_message("  u_L %12.6f;   u_R %12.6f\n", left.u[0], right.u[0]);
  log_message("  p_L %12.6f;   p_R %12.6f\n", left.p, right.p);
  float center = (pars.nx/2 + BC) * pars.dx;
  float wavevel;
  for (int i = BC; i<BC+pars.nx; i++){
    float x = (i+0.5)*pars.dx - center;
    float xovert = x / pars.tmax;
    riemann_solve(&left, &right, &grid[i].prim, xovert, &wavevel, /*dimension=*/0);
    grid[i].wavevel = wavevel;
  }


  /* write final */
  io_write_output(&outcount, step, t);



  printf("\n");
  printf("  Finished clean. Yay!\n");

  return(0);
}
