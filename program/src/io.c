/* IO routines */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cell.h"
#include "defines.h"
#include "gas.h" /* pstates */
#include "io.h"
#include "params.h"
#include "utils.h"

#if NDIM == 1
extern cell *grid;
#elif NDIM == 2
extern cell **grid;
#endif

extern params pars;

void io_read_cmdlineargs(int argc, char *argv[]) {
  /*--------------------------------------------------------
   * This function reads in the command line arguments and
   * stores them in the params struct
   *--------------------------------------------------------*/

  if (argc < 3) {
    throw_error("Too few arguments given. Run this program with ./hydro "
                "paramfile datafile\n");
  } else {
    strcpy(pars.paramfilename, argv[1]);
    strcpy(pars.datafilename, argv[2]);
  };
}

void io_read_ic_type(int *skip_lines) {
  /*--------------------------------------------------------
   * Start reading IC file, find what IC type we have
   *--------------------------------------------------------*/

  /* check whether file exists first */
  io_check_file_exists(pars.datafilename);

  /* open file */
  FILE *dat = fopen(pars.datafilename, "r");

  char varname[MAX_LINE_SIZE];
  char varvalue[MAX_LINE_SIZE];
  char tempbuff[MAX_LINE_SIZE];

  *skip_lines = 0;

  while (fgets(tempbuff, MAX_LINE_SIZE, dat)) {
    *skip_lines += 1;

    if (line_is_comment(tempbuff))
      continue;
    remove_trailing_comments(tempbuff);
    if (line_is_empty(tempbuff))
      continue;

    sscanf(tempbuff, "%20s = %56[^\n]\n", varname, varvalue);
    remove_whitespace(varname);
    remove_whitespace(varvalue);

    if (strcmp(varname, "filetype") == 0) {
      if (strcmp(varvalue, "two-state") == 0) {
        pars.twostate_ic = 1;
      } else if (strcmp(varvalue, "arbitrary") == 0) {
        pars.twostate_ic = 0;
      } else {
        throw_error(
            "while reading IC file type: I don't recognize file type \"%s\"\n",
            varvalue);
      }
      break;
    } else {
      throw_error(
          "while reading IC file type: Unrecongized line [error loc 1]\n    %s",
          tempbuff);
    }
  }

  /* If we have arbitrary IC file type, keep reading nx and ndim */
  char nx_read = 0;
  char ndim_read = 0;
  if (!pars.twostate_ic) {
    while (fgets(tempbuff, MAX_LINE_SIZE, dat)) {
      *skip_lines += 1;

      if (line_is_comment(tempbuff))
        continue;
      remove_trailing_comments(tempbuff);
      if (line_is_empty(tempbuff))
        continue;

      sscanf(tempbuff, "%20s = %56[^\n]\n", varname, varvalue);
      remove_whitespace(varname);
      remove_whitespace(varvalue);

      if (strcmp(varname, "nx") == 0) {
        pars.nx = atoi(varvalue);
        nx_read = 1;
      } else if (strcmp(varname, "ndim") == 0) {
        pars.ndim_ic = atoi(varvalue);
        ndim_read = 1;
      } else {
        throw_error("while reading IC file type: Unrecongized line [error loc "
                    "2]: \n    \"%s\" \n",
                    tempbuff);
      }

      if (nx_read && ndim_read)
        break;
    }
  }

  fclose(dat);

  debugmessage("In io_read_ic_type: got skiplines: %d", *skip_lines);
}

void io_read_ic_twostate(int skip) {
  /*--------------------------------------------------------
   * Read in initial conditions file, store read states.
   * This is for the two-state IC file format.
   * int skip: how many lines in the IC file to skip before
   * starting to read IC data
   *--------------------------------------------------------*/

  log_extra("Reading in IC");

  /* states to be read in */
  pstate left, right;

  gas_init_pstate(&left);
  gas_init_pstate(&right);

  /* bools to check we get complete data set */
  char rhol_set = 0;
  char ul_set = 0;
  char pl_set = 0;
  char rhor_set = 0;
  char ur_set = 0;
  char pr_set = 0;

  /* check whether file exists first */
  io_check_file_exists(pars.datafilename);

  /* open file */
  FILE *dat = fopen(pars.datafilename, "r");

  char varname[MAX_LINE_SIZE];
  char varvalue[MAX_LINE_SIZE];
  char tempbuff[MAX_LINE_SIZE];

  int i = 0;
  while (fgets(tempbuff, MAX_LINE_SIZE, dat)) {
    i += 1;
    /* debugmessage("i=%d, skip=%d, Got line: ||%s", i, skip, tempbuff); */
    if (i <= skip)
      continue; /* skip header */

    if (line_is_comment(tempbuff))
      continue;
    remove_trailing_comments(tempbuff);
    if (line_is_empty(tempbuff))
      continue;
    if (!check_name_equal_value_present(tempbuff))
      continue;

    sscanf(tempbuff, "%20s = %56[^\n]\n", varname, varvalue);
    remove_whitespace(varname);
    remove_whitespace(varvalue);

    if (strcmp(varname, "rho_L") == 0) {
      left.rho = atof(varvalue);
      rhol_set = 1;
    } else if (strcmp(varname, "u_L") == 0) {
      left.u[0] = atof(varvalue);
      ul_set = 1;
    } else if (strcmp(varname, "p_L") == 0) {
      left.p = atof(varvalue);
      pl_set = 1;
    } else if (strcmp(varname, "rho_R") == 0) {
      right.rho = atof(varvalue);
      rhor_set = 1;
    } else if (strcmp(varname, "u_R") == 0) {
      right.u[0] = atof(varvalue);
      ur_set = 1;
    } else if (strcmp(varname, "p_R") == 0) {
      right.p = atof(varvalue);
      pr_set = 1;
    } else {
      log_message("ATTENTION: Unrecongized data : \"%s\" | \"%s\"\n", varname,
                  varvalue);
    }
  }

  fclose(dat);

  if (!rhol_set)
    throw_error("rho left is not given in IC file");
  if (!rhor_set)
    throw_error("rho right is not given in IC file");
  if (!ul_set)
    throw_error("u left is not given in IC file");
  if (!ur_set)
    throw_error("u right is not given in IC file");
  if (!pl_set)
    throw_error("u left is not given in IC file");
  if (!pr_set)
    throw_error("u right is not given in IC file");

    /* Now write the data in the actual grid */

#if NDIM == 1
  for (int i = BC; i < pars.nx / 2 + BC; i++) {
    grid[i].prim.rho = left.rho;
    grid[i].prim.u[0] = left.u[0];
    grid[i].prim.p = left.p;
  }
  for (int i = pars.nx / 2 + BC; i < pars.nx + BC; i++) {
    grid[i].prim.rho = right.rho;
    grid[i].prim.u[0] = right.u[0];
    grid[i].prim.p = right.p;
  }

#elif NDIM == 2

  for (int j = BC; j < pars.nx + BC; j++) {

    for (int i = BC; i < pars.nx / 2 + BC; i++) {
      grid[i][j].prim.rho = left.rho;
      grid[i][j].prim.u[0] = left.u[0];
      grid[i][j].prim.u[1] = 0;
      grid[i][j].prim.p = left.p;
    }
    for (int i = pars.nx / 2 + BC; i < pars.nx + BC; i++) {
      grid[i][j].prim.rho = right.rho;
      grid[i][j].prim.u[0] = right.u[0];
      grid[i][j].prim.u[1] = 0;
      grid[i][j].prim.p = right.p;
    }
  }

#endif

  log_message("The initial discontinuity is at x = %8.5lf\n",
              (pars.nx / 2) * pars.dx);
}

void io_read_ic_arbitrary(int skip) {
  /*--------------------------------------------------------
   * Read in initial conditions file, store read states.
   * This is for the arbitrary IC file format.
   * int skip: how many lines in the IC file to skip before
   * starting to read IC data
   *--------------------------------------------------------*/

  log_extra("Reading in IC");

  /* check whether file exists first */
  io_check_file_exists(pars.datafilename);

  /* open file */
  FILE *dat = fopen(pars.datafilename, "r");

  char tempbuff[MAX_LINE_SIZE];

  int i = BC;
#if NDIM > 1
  int j = BC;
#endif
  int sc = 0;
  int counter = 0;
  while (fgets(tempbuff, MAX_LINE_SIZE, dat)) {
    sc += 1;
    if (sc <= skip)
      continue; /* skip header */

    if (line_is_comment(tempbuff))
      continue;
    remove_trailing_comments(tempbuff);
    if (line_is_empty(tempbuff))
      continue;

#if NDIM == 1
    float rho, u, p;
    check_number_of_columns_IC(tempbuff, 3);
    sscanf(tempbuff, "%f %f %f\n", &rho, &u, &p);

    grid[i].prim.rho = rho;
    grid[i].prim.u[0] = u;
    grid[i].prim.u[1] = 0;
    grid[i].prim.p = p;

    i += 1;

#elif NDIM == 2
    float rho, ux, uy, p;
    check_number_of_columns_IC(tempbuff, 4);
    sscanf(tempbuff, "%f %f %f %f\n", &rho, &ux, &uy, &p);

    grid[i][j].prim.rho = rho;
    grid[i][j].prim.u[0] = ux;
    grid[i][j].prim.u[1] = uy;
    grid[i][j].prim.p = p;

    i += 1;
    if (i == pars.nx + BC) {
      i = BC;
      j += 1;
    }

#endif
    counter += 1;

    /* safety measure */
#if NDIM == 1
    if (counter == pars.nx)
      break;
#elif NDIM == 2
    if (counter == pars.nx * pars.nx)
      break;
#endif
  }

  /* another safety mesure */
#if NDIM == 1
  if (counter != pars.nx) {
    throw_error("In io_read_ic_arbitrary:I didn't get the proper number of IC "
                "cell data.\n"
                "Expected nx = %d, got %d",
                pars.nx, counter);
  }
#elif NDIM == 2
  if (counter != pars.nx * pars.nx) {
    throw_error("In io_read_ic_arbitrary: I didn't get the proper number of IC "
                "cell data.\n"
                "Expected nx^2 = %d, got %d",
                pars.nx * pars.nx, counter);
  }
#endif
  else {
    while (fgets(tempbuff, MAX_LINE_SIZE, dat)) {
      if (line_is_comment(tempbuff))
        continue;
      remove_trailing_comments(tempbuff);
      if (line_is_empty(tempbuff))
        continue;

      /* if you arrived at this point, you have something that is not a comment
       * nor empty line even though we have all the data that we need */
      throw_error("In io_read_ic_arbitrary: I read enough data according to "
                  "given nx, but there is still more stuff in the IC file.\n"
                  "The line I read was: %s",
                  tempbuff);
    }
  }

  fclose(dat);
}

void io_read_paramfile() {
  /*------------------------------------------------------------*/
  /* Read in parameter file, store read in global parameters.   */
  /*------------------------------------------------------------*/

  /* check whether file exists first */
  io_check_file_exists(pars.paramfilename);

  /* open file */
  FILE *par = fopen(pars.paramfilename, "r");

  char varname[MAX_LINE_SIZE];
  char varvalue[MAX_LINE_SIZE];
  char tempbuff[MAX_LINE_SIZE];

  while (fgets(tempbuff, MAX_LINE_SIZE, par)) {

    if (line_is_comment(tempbuff))
      continue;
    remove_trailing_comments(tempbuff);
    if (line_is_empty(tempbuff))
      continue;
    if (!check_name_equal_value_present(tempbuff))
      continue;

    sscanf(tempbuff, "%20s = %56[^\n]\n", varname, varvalue);
    remove_whitespace(varname);
    remove_whitespace(varvalue);

    if (strcmp(varname, "verbose") == 0) {
      pars.verbose = atoi(varvalue);
    } else if (strcmp(varname, "nstep_log") == 0) {
      pars.nstep_log = atoi(varvalue);
    } else if (strcmp(varname, "nsteps") == 0) {
      pars.nsteps = atoi(varvalue);
    } else if (strcmp(varname, "tmax") == 0) {
      pars.tmax = atof(varvalue);
    } else if (strcmp(varname, "nx") == 0) {
      pars.nx = atoi(varvalue);
    } else if (strcmp(varname, "ccfl") == 0) {
      pars.ccfl = atof(varvalue);
    } else if (strcmp(varname, "force_dt") == 0) {
      pars.force_dt = atof(varvalue);
    } else if (strcmp(varname, "boundary") == 0) {
      pars.boundary = atoi(varvalue);
    } else if (strcmp(varname, "foutput") == 0) {
      pars.foutput = atoi(varvalue);
    } else if (strcmp(varname, "dt_out") == 0) {
      pars.dt_out = atof(varvalue);
    } else if (strcmp(varname, "basename") == 0) {
      if (!line_is_empty(varvalue)) {
        strcpy(pars.outputfilename, varvalue);
      }
    } else if (strcmp(varname, "toutfile") == 0) {
      if (!line_is_empty(varvalue)) {
        strcpy(pars.toutfilename, varvalue);
        pars.use_toutfile = 1;
      }
    } else if (strcmp(varname, "src_const_acc_x") == 0) {
      pars.src_const_acc_x = atof(varvalue);
      pars.sources_are_read = 1;
    } else if (strcmp(varname, "src_const_acc_y") == 0) {
      pars.src_const_acc_y = atof(varvalue);
      pars.sources_are_read = 1;
    } else if (strcmp(varname, "src_const_acc_r") == 0) {
      pars.src_const_acc_r = atof(varvalue);
      pars.sources_are_read = 1;
    } else {
      log_message("ATTENTION: Unrecongized parameter : \"%s\"\n", varname);
    }
  }

  fclose(par);
}

void io_read_toutfile() {
  /*------------------------------------------------------------*/
  /* Read in parameter file, store read in global parameters.   */
  /*------------------------------------------------------------*/

  /* check whether file exists first */
  io_check_file_exists(pars.toutfilename);

  /* open file */
  FILE *par = fopen(pars.toutfilename, "r");

  char tempbuff[MAX_LINE_SIZE];

  int nlines = 0;
  float past_value = 0;
  float value = 0;

  /* get how many lines we have */
  while (fgets(tempbuff, MAX_LINE_SIZE, par)) {

    if (line_is_comment(tempbuff))
      continue;
    remove_trailing_comments(tempbuff);
    if (line_is_empty(tempbuff))
      continue;

    nlines += 1;
    sscanf(tempbuff, "%f\n", &value);

    if (nlines > 1) {
      if (value <= past_value)
        throw_error("While reading output time file '%s':"
                    " Got output time %f <= previous output time %f",
                    pars.toutfilename, value, past_value);
    }
    past_value = value;
  }
  fclose(par);

  /* Now allocate the array and read in the stuff */

  pars.outputtimes = malloc(nlines * sizeof(float));
  pars.noutput_tot = nlines;
  nlines = 0;

  par = fopen(pars.toutfilename, "r");
  while (fgets(tempbuff, MAX_LINE_SIZE, par)) {

    if (line_is_comment(tempbuff))
      continue;
    remove_trailing_comments(tempbuff);
    if (line_is_empty(tempbuff))
      continue;

    sscanf(tempbuff, "%f\n", &pars.outputtimes[nlines]);
    nlines += 1;
  }
  fclose(par);
}

void io_write_output(int *outstep, int step, float t) {
  /*----------------------------------------*/
  /* Write output of step at time t.
   *
   * outstep: Current output number
   * step:  Current step of the simulation
   * t:     Current time of the simulation
   *----------------------------------------*/

  if (*outstep > 9999)
    throw_error("I'm not made to write outputs > 9999\n");

  /* generate output filename for this snapshot */

  char filename[MAX_FNAME_SIZE] = "";
  char snapnrstr[5] = "";

  strcpy(filename, pars.outputfilename);
  sprintf(snapnrstr, "%04d", *outstep);
  strcat(filename, "-");
  strcat(filename, snapnrstr);
  strcat(filename, ".out");

  log_extra("Dumping output to %s for t= %g", filename, t);

  /*-------------------*/
  /* Write output file */
  /*-------------------*/

  FILE *outfilep = fopen(filename, "w");
  fprintf(outfilep, "# ndim = %2d\n", NDIM);
  fprintf(outfilep, "# nx = %10d\n", pars.nx);
  fprintf(outfilep, "# t = %12.6lf\n", t);
  fprintf(outfilep, "# nsteps = %12d\n", step);

#if NDIM == 1
  fprintf(outfilep, "# %12s %12s %12s %12s\n", "x", "rho", "u", "p");
  for (int i = BC; i < pars.nx + BC; i++) {
    cell c = grid[i];
    pstate s = c.prim;
    fprintf(outfilep, "%12.6e %12.6e %12.6e %12.6e\n", c.x, s.rho, s.u[0], s.p);
  }

#elif NDIM == 2

  fprintf(outfilep, "# %12s %12s %12s %12s %12s %12s\n", "x", "y", "rho", "u_x",
          "u_y", "p");
  for (int j = BC; j < pars.nx + BC; j++) {
    for (int i = BC; i < pars.nx + BC; i++) {
      cell c = grid[i][j];
      pstate s = c.prim;
      fprintf(outfilep, "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n", c.x, c.y,
              s.rho, s.u[0], s.u[1], s.p);
    }
  }

#endif
  fclose(outfilep);

  /* raise output step number */
  *outstep += 1;
}

int io_is_output_step(float t, float *dt, int step) {
  /* -------------------------------------------------------------------------------------------------
   * Check whether we should be writing an output in this time step. Returns 1
   * if true, 0 otherwise. If necessary, reduces size of dt so that it fits
   * required output time exactly.
   *
   * t:     current time of sim
   * dt:    current time step size of sim
   * step:  current step of sim
   * -------------------------------------------------------------------------------------------------
   */

  debugmessage(
      "Checking whether we need to limit the timestep for output. t=%g, dt=%g",
      t, *dt);

  if (pars.use_toutfile) { /* if we have a toutfile or dt_out in params is given
                            */
    if (pars.noutput >= pars.noutput_tot)
      return (0); /* final output will be dumped anyhow */
    float tnext = pars.outputtimes[pars.noutput];
    if (t + *dt >= tnext) {
      debugmessage("Overwriting dt from %f to %f such that t+dt=%f", *dt,
                   tnext - t, tnext);
      *dt = tnext - t;
      pars.noutput += 1;
      return (1);
    }
  }

  /* do timestep limiting first! */
  /* this case can happen if dt_out = 0 and foutput = 0! */
  if (pars.tmax > 0 && pars.tmax < t + *dt) {
    debugmessage("Overwriting dt from %f to %f such that t+dt=%f", *dt,
                 pars.tmax - t, pars.tmax);
    *dt = pars.tmax - t;
    return (
        0); /* write up as not output step because last output will be dumped */
  }

  if (pars.foutput > 0) {
    if (step % pars.foutput == 0) {
      return (1);
    } else {
      return (0);
    }
  }

  return (0);
}

void io_check_file_exists(char *fname) {
  /* -------------------------------------------------------- */
  /* Check whether a file exists. If it doesn't, exit.        */
  /* -------------------------------------------------------- */

  FILE *f = fopen(fname, "r");

  /* check if file exists */
  if (f == NULL) {
    throw_error("file '%s' not found.\n", fname);
  } else {
    fclose(f);
  }
}

int line_is_empty(char *line) {
  /* --------------------------------- */
  /* Check whether this line is empty, */
  /* i.e. only whitespaces or newlines.*/
  /* returns 1 if true, 0 otherwise.   */
  /* assumes line is MAX_LINE_SIZE     */
  /* --------------------------------- */

  int isempty = 0;

  for (int i = 0; i < MAX_LINE_SIZE; i++) {
    if (line[i] != ' ') {
      if (line[i] == '\n')
        isempty = 1;
      break;
    }
  }
  return (isempty);
}

int line_is_comment(char *line) {
  /* --------------------------------------
   * Check whether the given line string is
   * a comment, i.e. starts with // or
   * <slash>*
   * -------------------------------------- */

  /* initialize firsttwo explicily for valgrind */
  char firsttwo[3] = {0, 0, 0};

  strncpy(firsttwo, line, 2);

  /* strcmp returns 0 if strings are equal */
  if (!strcmp(firsttwo, "//") || !strcmp(firsttwo, "/*")) {
    return (1);
  } else {
    return (0);
  }
}

void remove_whitespace(char *line) {
  /*---------------------------------------------------------
   * remove heading and trailing whitespaces
   * --------------------------------------------------------*/

  int start = 0;
  int stop = strlen(line);

  /* find first non-whitespace character */
  for (int i = 0; i < MAX_LINE_SIZE; i++) {
    if ((line[i] != ' ') && (line[i] != '\t')) {
      start = i;
      break;
    }
  }

  /* find last non-whitespace character */
  for (int i = 0; i < stop; i++) {
    if ((line[stop - i - 1] != ' ') && (line[stop - i - 1] != '\t')) {
      stop = stop - i - 1;
      break;
    }
  }

  char newline[MAX_LINE_SIZE];
  strncpy(newline, line + start, stop - start + 1);
  newline[stop - start + 1] = '\0';
  strcpy(line, newline);
}

void remove_trailing_comments(char *line) {
  /*---------------------------------------------------------
   * Check whether there are trailing comments in this line
   * and if so, remove them.
   * --------------------------------------------------------*/

  for (int i = 0; i < MAX_LINE_SIZE - 2; i++) {
    /* -2: 1 for \0 char, 1 because comment is always 2 characters long,
     * and I only check for the first.*/
    if (line[i] == '\0') {
      break;
    } else if (line[i] == '/') {
      char twochars[3];
      strncpy(twochars, line + i, 2);
      if (line_is_comment(twochars)) {
        line[i] = '\n';
        line[i + 1] = '\0';
        break;
      }
    }
  }
}

int check_name_equal_value_present(char *line) {
  /* ----------------------------------------------------------------
   * Check that the line you're reading has the correct number of
   * columns, delimited by an equality sign.
   * We expect <name> = <param>
   *
   * returns 1 if that is the case, 0 if not.
   * ----------------------------------------------------------------*/

  int pos = 0;
  int check = 0;
  int len = strlen(line);

  /* first check that we have an equal sign in there */
  for (int i = 0; i < len; i++) {
    if (line[i] == '=') {
      pos = i;
      break;
    }
  }
  if (pos == 0)
    return (0); /* '=' is either in first place or not present. */

  /* now check that we have non-whitespace characters */
  for (int i = 0; i < pos; i++) {
    if (line[i] != ' ') {
      check = 1;
      break;
    }
  }

  if (!check)
    return (0);
  check = 0;

  /* now check that we have non-whitespace characters after the equal sign */
  for (int i = pos + 1; i < len; i++) {
    if (line[i] != ' ') {
      check = 1;
      break;
    }
  }

  if (!check)
    return (0);
  return (1);
}

void check_number_of_columns_IC(char *line, int should) {
  /* ----------------------------------------------------------------
   * Check that the line you're reading has the correct number of
   * columns, delimited by an empty space. int should defines the
   * expected number of columns.
   * ----------------------------------------------------------------*/

  int splits = 0;
  char cpy[MAX_LINE_SIZE];
  strcpy(cpy, line);

  /* remove trailing whitespace */
  int i = strlen(cpy);
  while (cpy[i - 1] == ' ' || cpy[i - 1] == '\n') {
    cpy[i - 1] = '\0';
    i -= 1;
  }

  char *split_token = strtok(cpy, " ");

  if (split_token != NULL) {
    while (split_token != NULL) {
      splits += 1;
      split_token = strtok(NULL, " ");
    }
  }
  if (splits != should) {
    throw_error("I expected %d columns, I got %d; Line was %s", should, splits,
                line);
  }
}
