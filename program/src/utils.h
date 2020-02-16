/* Misc utils that don't fit anywhere else*/

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef UTILS_H
#define UTILS_H

#include "defines.h"

/* stuff that doesn't fit anywhere else */
void print_header();
void print_compile_defines();
void utils_get_macro_strings(char* solver, char* riemann, char *limiter);


/* helper functions */
void log_message(const char *format, ...);
void debugmessage(const char *format, ...);
void log_extra(const char *format, ...);
void throw_error(const char *format, ...);
void printbool(int boolean);

#endif
