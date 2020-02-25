/* Cell related stuff */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef CELL_H
#define CELL_H


#include "gas.h"

typedef struct{

  int id;

  float x;     /* x position of cell center */
  float y;     /* y position of cell center */

  pstate prim;  /* primitive state vector */
  cstate cons;  /* conserved state vector */

  pstate flux;  /* fluxes of primitive variables */

  float wavevel; /* highest wave velocity */

} cell;



#if NDIM == 1
extern cell* grid;
#elif NDIM == 2
extern cell** grid;
#endif


void cell_init_cell(cell* c);
void cell_init_grid();
void cell_set_boundary();
void cell_real_to_ghost(cell** realL, cell** realR, cell** ghostL, cell** ghostR);
void cell_copy_boundary_data(cell* real, cell* ghost);
void cell_copy_boundary_data_reflective(cell* real, cell* ghost);
void cell_reset_fluxes();

float cell_get_total_mass();

void cell_print_grid(char field[4]);
void cell_print_grid_part(char field[4], int* limits);
void cell_get_ij(cell* c, int* i, int* j);
char* cell_get_index_string(cell* c);

#endif
