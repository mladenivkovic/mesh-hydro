/* Cell and grid related stuff */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef CELL_H
#define CELL_H

#include "defines.h"
#include "gas.h"

typedef struct {

  int id;

  float x; /* x position of cell center */
  float y; /* y position of cell center */

  pstate prim; /* primitive state vector */
  cstate cons; /* conserved state vector */

  pstate pflux; /* fluxes of primitive variables */
  cstate cflux; /* fluxes of conserved variables */

  float acc[2]; /* acceleration from external sources */

#if SOLVER == WAF
  /* The following is needed for the WAF scheme */
  float Sk[3];              /* wave speeds of the Riemann solution */
  cstate riemann_fluxes[4]; /* the emerging fluxes: F_L, F*_L, F*_R, F_R */
  float delta_q[3];         /* difference in densities between each wave */
#endif
#if SOLVER == MUSCL
  /* The following is needed for the MUSCL scheme */
  cstate ULmid; /* updated left extrapolated boundary value */
  cstate URmid; /* updated right extrapolated boundary value */
#endif

} cell;

#if NDIM == 1
extern cell *grid;
#elif NDIM == 2
extern cell **grid;
#endif

void cell_init_cell(cell *c);
void cell_init_grid();
void cell_set_boundary();
void cell_real_to_ghost(cell **realL, cell **realR, cell **ghostL,
                        cell **ghostR, int dimension);
void cell_copy_boundary_data(cell *real, cell *ghost);
void cell_copy_boundary_data_reflective(cell *real, cell *ghost, int dimension);
void cell_reset_fluxes();
void cell_get_pstates_from_cstates();
void cell_get_cstates_from_pstates();

float cell_get_total_mass();

void cell_print_grid(char field[4]);
void cell_print_grid_part(char field[4], int *limits);
void cell_get_ij(cell *c, int *i, int *j);
char *cell_get_index_string(cell *c);

#endif
