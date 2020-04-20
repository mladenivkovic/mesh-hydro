/* IO routines */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */


#ifndef IO_H
#define IO_H

#include "defines.h"

void io_read_cmdlineargs(int argc, char* argv[]);
void io_read_ic_type(int* skip_lines);
void io_read_ic_twostate();
void io_read_ic_arbitrary(int skip);
void io_read_paramfile();
void io_read_toutfile();
void io_write_output(int *outstep, int step, float t);
int io_is_output_step(float t, float* dt, int step);



void io_check_file_exists(char* fname);
int line_is_empty(char* line);
int line_is_comment(char* line);
void remove_whitespace(char* line);
void remove_trailing_comments(char* line);
int check_name_equal_value_present(char* line);
void check_number_of_columns_IC(char* line, int should);

#endif
