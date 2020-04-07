/* Misc utils that don't fit anywhere else*/

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"
#include "utils.h"

extern params pars;




void print_header(){
  /*----------------------*/
  /* Print a nice header  */
  /*----------------------*/
  printf("===============================================================================================\n");
  printf(" ,--.   ,--.,------. ,---.  ,--.  ,--.    ,--.  ,--.,--.   ,--.,------.  ,------.  ,-----.  \n"); 
  printf(" |   `.'   ||  .---''   .-' |  '--'  |    |  '--'  | \\  `.'  / |  .-.  \\ |  .--. ''  .-.  ' \n");
  printf(" |  |'.'|  ||  `--, `.  `-. |  .--.  |    |  .--.  |  '.    /  |  |  \\  :|  '--'.'|  | |  | \n");
  printf(" |  |   |  ||  `---..-'    ||  |  |  |    |  |  |  |    |  |   |  '--'  /|  |\\  \\ '  '-'  ' \n");
  printf(" `--'   `--'`------'`-----' `--'  `--'    `--'  `--'    `--'   `-------' `--' '--' `-----'  \n");
  printf("===============================================================================================\n");
}





void print_compile_defines(){
  /* ----------------------------------- */
  /* print compile time definitions      */
  /* ----------------------------------- */


  char solver[80];
  char riemann[80];
  char limiter[80];

  utils_get_macro_strings(solver, riemann, limiter);

  log_message("-----------------------------------------------------------------------------------------\n");
  log_message("\n");
  log_message("Compile time parameters are:\n");
  log_message("\n");

  log_message("Compile time:                "STR(COMPDATE)"\n");
  log_message("Dimensions:                  "STR(NDIM)"\n");
  log_message("Hydro solver:                %s\n", solver);
  log_message("Riemann solver:              %s\n", riemann);
  log_message("Limiter:                     %s\n", limiter);
}





void utils_get_macro_strings(char* solver, char* riemann, char* limiter){
  /* ------------------------------------------------------------------------ 
   * Get string names for the solver, riemann solver, and limiter in use.
   * ------------------------------------------------------------------------ */

#if SOLVER == ADVECTION_PWCONST
  strcpy(solver, "ADVECTION_PWCONST");
#elif SOLVER == ADVECTION_PWLIN
  strcpy(solver, "ADVECTION_PWLIN");
#elif SOLVER == GODUNOV
  strcpy(solver, "GODUNOV");
#elif SOLVER == ADVECTION_WAF
  strcpy(solver, "ADVECTION_WAF");
#elif SOLVER == WAF
  strcpy(solver, "WAF");
#elif SOLVER == MUSCL
  strcpy(solver, "MUSCL");
#elif SOLVER == NONE
  strcpy(solver, "");
#endif


#if RIEMANN == NONE
  strcpy(riemann, "NONE");
#elif RIEMANN == EXACT
  strcpy(riemann, "EXACT");
#elif RIEMANN == HLLC
  strcpy(riemann, "HLLC");
#elif RIEMANN == HLL
  strcpy(riemann, "HLL");
#elif RIEMANN == TRRS
  strcpy(riemann, "TRRS");
#elif RIEMANN == TSRS
  strcpy(riemann, "TSRS");
#else
  strcpy(riemann, "unknown");
#endif


#if LIMITER == NONE
  strcpy(limiter, "NONE");
#elif LIMITER == MINMOD
  strcpy(limiter, "MINMOD");
#elif LIMITER == SUPERBEE
  strcpy(limiter, "SUPERBEE");
#elif LIMITER == VANLEER
  strcpy(limiter, "VAN_LEER");
#elif LIMITER == MC
  strcpy(limiter, "MC");
#endif
}




void log_message(const char *format, ...){
  /*------------------------------------------------------
   * just prepends [LOG] to printing. Use like printf()
   *------------------------------------------------------*/

  if (pars.verbose < 1) return;

  printf("%-12s", "[LOG] ");

  va_list arg; /* from stdarg.h; va_list type defined in stdarg.h */
  int done;

  va_start(arg, format); /* initialises arg variable, starting with "format" */
                         /* variable given as an argument to print() */

  done = vfprintf(stdout, format, arg); /* call the formatting and printing */
                                        /* function that printf also uses */

  va_end(arg); /* do whatever cleanup is necessary */

  if (done < 0) throw_error("My own log_message() function exited with error code %d", done);
}



void debugmessage(const char *format, ...){
  /*---------------------------------------------------------------
   * if verbose is 3 or higher, write whatever you want to write 
   * to screen. Use it like you use printf(), except this function 
   * will add a newline by itself :)
   *---------------------------------------------------------------*/

  if (pars.verbose < 3) return;

  printf("%-12s", "[DEBUGGING] ");

  va_list arg; /* from stdarg.h; va_list type defined in stdarg.h */
  int done;

  va_start(arg, format); /* initialises arg variable, starting with "format" */
                         /* variable given as an argument to print() */

  done = vfprintf(stdout, format, arg); /* call the formatting and printing */
                                        /* function that printf also uses */

  va_end(arg); /* do whatever cleanup is necessary */

  printf("\n");

  if (done < 0) throw_error("My own debugmessage() function exited with error code %d", done);
}



void log_extra(const char *format, ...){
  /*-----------------------------------------------
   * if verbose is 2 or higher, write whatever you
   * want to write to screen. Use it like you use
   * printf(), except this function will add a 
   * newline by itself :) .
   *-----------------------------------------------*/

  if (pars.verbose < 2) return;

  printf("%-12s", "[EXTRA] ");

  va_list arg; /* from stdarg.h; va_list type defined in stdarg.h */
  int done;

  va_start(arg, format); /* initialises arg variable, starting with "format" */
                         /* variable given as an argument to print() */

  done = vfprintf(stdout, format, arg); /* call the formatting and printing */
                                        /* function that printf also uses */

  va_end(arg); /* do whatever cleanup is necessary */

  printf("\n");

  if (done < 0) throw_error("My own log_extra() function exited with error code %d", done);
}



void throw_error(const char *format, ...){
  /*-----------------------------------------------
   * Print a formatted error message to screen.
   * Use it like you use printf(), except this 
   * function will add a newline by itself :) .
   * Then it will exit.
   *-----------------------------------------------*/

  printf("%-12s", "[ERROR]");

  va_list arg; /* from stdarg.h; va_list type defined in stdarg.h */
  int done;

  va_start(arg, format); /* initialises arg variable, starting with "format" */
                         /* variable given as an argument to print() */

  done = vfprintf(stdout, format, arg); /* call the formatting and printing */
                                        /* function that printf also uses */

  va_end(arg); /* do whatever cleanup is necessary */

  printf("\n"); /* always end with a newline! :) */

  if (done < 0) throw_error("ERROR: your errormessage() function exited with error code %d\n", done);

  printf("GG yall, I'm out\n");
  exit(1);
}



void printbool(int boolean){
  /*-------------------------------------*/
  /* prints "True" or "False"            */
  /*-------------------------------------*/

  if (boolean){
    printf("true");
  } else {
    printf("false");
  }
}
