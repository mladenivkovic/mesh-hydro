/* =======================================================
 * Use the code to only solve a Riemann problem.
 * run with ./riemann parameterfile.txt ICfile.dat
 * See the README.md in parent directory on how the
 * parameter and IC file can look like.
 *
 * Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com
 * ======================================================= */

#include <stdio.h>
#include <stdlib.h>

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
int main(int argc, char *argv[]) {
  /* ====================================== */

  /* initialize parameters struct */
  params_init_defaults();

  /* read cmdline args and parameter file */
  io_read_cmdlineargs(argc, argv);
  io_read_paramfile();

  /* check which IC type we have */
  int skiplines_ic =
      0; /* how many lines to skip next time you open the IC file */
  io_read_ic_type(&skiplines_ic);

  if (!pars.twostate_ic) {
    printf("I can only solve twostate-ICs.");
    exit(2);
  }

  params_generate_riemann_output_filename(); /* create a different file name
                                                than the default */
  params_init_derived();  /* process the parameters you got. */
  params_check_riemann(); /* check whether we can work with this setup. */

  /* print / announce stuff for logging */
  print_compile_defines();
  params_print_log();

  /* initialize the grid */
  cell_init_grid();

  /* read in the full IC file */
  io_read_ic_twostate(skiplines_ic);

  /* get the left and right state of the original problem */
  /* copy them, otherwise you will overwrite them! */
  pstate left;
  gas_init_pstate(&left);
  left.rho = grid[pars.nx / 2 + BC - 1].prim.rho;
  left.u[0] = grid[pars.nx / 2 + BC - 1].prim.u[0];
  left.u[1] = grid[pars.nx / 2 + BC - 1].prim.u[1];
  left.p = grid[pars.nx / 2 + BC - 1].prim.p;
  pstate right;
  gas_init_pstate(&right);
  right.rho = grid[pars.nx / 2 + BC].prim.rho;
  right.u[0] = grid[pars.nx / 2 + BC].prim.u[0];
  right.u[1] = grid[pars.nx / 2 + BC].prim.u[1];
  right.p = grid[pars.nx / 2 + BC].prim.p;

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

  /* center: where the initial separation of the states is.*/
  float center = (pars.nx / 2) * pars.dx;
  for (int i = BC; i < BC + pars.nx; i++) {
    float x = (i - BC + 0.5) * pars.dx - center;
    float xovert = x / pars.tmax;
#if RIEMANN == HLLC
    riemann_solve_hllc_state(&left, &right, &grid[i].prim, xovert,
                             /*dimension=*/0);
#else
    riemann_solve(&left, &right, &grid[i].prim, xovert, /*dimension=*/0);
#endif
  }

  /* write final */
  io_write_output(&outcount, step, t);

  printf("\n");
  printf("  Finished clean. Yay!\n");

  return (0);
}
