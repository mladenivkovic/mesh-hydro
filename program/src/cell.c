/* Cell related stuff */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "cell.h"
#include "defines.h"
#include "gas.h"
#include "params.h"
#include "utils.h"




#if NDIM == 1
extern cell *grid;
#elif NDIM == 2
extern cell **grid;
#endif

extern params pars;





void cell_init_cell(cell * c){
  /*---------------------------------------*/
  /* Initialize/reset the values of a cell */
  /*---------------------------------------*/

  c->id = 0;
  c->x = 0;
  c->y = 0;
  
  gas_init_pstate(&(c->prim));
  gas_init_pstate(&(c->flux));
  gas_init_cstate(&(c->cons));
}





void cell_init_grid(){
  /*----------------------------------------------------
   * Initialize the grid: Allocate memory for the cells,
   * initialize them, then distribute their position.
   * Convention: 
   *    grid[0][0]       is lower left corner
   *    grid[nxtot-1][0] is lower right corner
   *    grid[0][nxtot-1] is upper left corner
   *----------------------------------------------------*/

  log_extra("Initializing grid; ndim=%d, nx=%d", NDIM, pars.nx);

#if NDIM == 1
  
  grid = malloc(pars.nxtot * sizeof(cell));
  
  for (int i = 0; i < pars.nxtot; i++){
    cell_init_cell(&grid[i]);
    grid[i].x = (i-BC+0.5)*pars.dx;
    grid[i].id = i;
  }

#elif NDIM == 2
  
  grid = malloc(pars.nxtot * sizeof(cell*));
  
  for (int i = 0; i < pars.nxtot; i++){
    grid[i] = malloc(pars.nxtot * sizeof(cell));
    for (int j = 0; j < pars.nxtot; j++){
      cell_init_cell(&grid[i][j]);
    }
  }

  for (int j = 0; j < pars.nxtot; j++){
    for (int i = 0; i < pars.nxtot; i++){
      grid[i][j].x = (i-BC+0.5)*pars.dx;
      grid[i][j].y = (j-BC+0.5)*pars.dx;
      grid[i][j].id = i + j*pars.nxtot;
    }
  }
#endif
}




void cell_set_boundary(){
  /* -------------------------------------------------------
   * enforce boundary conditions. 
   * This function only picks out the pairs of real
   * and ghost cells in a row or column and then 
   * calls the function that actually copies the data.
   * ------------------------------------------------------- */

  debugmessage("Setting boundary conditions");

  cell* realL[BC];
  cell* realR[BC];
  cell* ghostL[BC];
  cell* ghostR[BC];

#if NDIM == 1

  for(int i = 0; i < BC; i++){
    realL[i] = &(grid[BC+i]);
    realR[i] = &(grid[pars.nx+i]); /* = last index of a real cell  - BC + (i + 1) */
    ghostL[i] = &(grid[i]);
    ghostR[i] = &(grid[pars.nx+BC+i]);
  }
  cell_real_to_ghost(realL, realR, ghostL, ghostR);

#elif NDIM == 2

  /* first do all left-right boundaries */
  for (int j = 0; j < pars.nx + BCTOT; j++){
    for(int i = 0; i < BC; i++){
      realL[i] = &(grid[BC+i][j]);
      realR[i] = &(grid[pars.nx+i][j]); /* = last index of a real cell  - BC + (i + 1) */
      ghostL[i] = &(grid[i][j]);
      ghostR[i] = &(grid[pars.nx+BC+i][j]);
    }
    cell_real_to_ghost(realL, realR, ghostL, ghostR);
  }

  /* now do all upper-lower boundaries */
  for (int i = 0; i < pars.nx + BCTOT; i++){
    for(int j = 0; j < BC; j++){
      realL[j] = &(grid[i][BC+j]);
      realR[j] = &(grid[i][pars.nx+j]); /* = last index of a real cell  - BC + (j + 1) */
      ghostL[j] = &(grid[i][j]);
      ghostR[j] = &(grid[i][pars.nx+BC+j]);
    }
    cell_real_to_ghost(realL, realR, ghostL, ghostR);
  }

#endif
}





void cell_real_to_ghost(cell** realL, cell** realR, cell** ghostL, cell** ghostR){
  /* ------------------------------------------------------------
   * apply the boundary conditions from real to ghost cells
   *
   * realL: array of pointers to real cells with lowest index
   * realR: array of pointers to real cells with highest index
   * ghostL: array of pointers to ghost cells with lowest index
   * ghostR: array of pointers to ghost cells with highest index
   *
   * all arguments are arrays of size BC, defined in defines.h
   * lowest array index is also lowest index of cell in grid
   * ------------------------------------------------------------ */

  if (pars.boundary == 0){
    /* periodic boundary conditions */
    /* ---------------------------- */
    for (int i = 0; i < BC; i++){
      cell_copy_boundary_data(realL[i], ghostR[i]);
      cell_copy_boundary_data(realR[i], ghostL[i]);
    }


  } else if (pars.boundary == 1){
    /* reflective boundary conditions */
    /* ------------------------------ */
    for (int i = 0; i < BC; i++){
      cell_copy_boundary_data_reflective(realL[i], ghostL[BC-1-i]);
      cell_copy_boundary_data_reflective(realR[i], ghostR[BC-1-i]);
    }


  } else if (pars.boundary == 2){
    /* transmissive boundary conditions */
    /* -------------------------------- */
    for (int i = 0; i < BC; i++){
      cell_copy_boundary_data(realL[0], ghostL[i]);
      cell_copy_boundary_data(realR[BC-1], ghostR[i]);
    }

  }
}







void cell_copy_boundary_data(cell* real, cell* ghost){
  /* ---------------------------------------------------------
   * Copies the actual data needed for boundaries from a real 
   * cell to a ghost cell
   * --------------------------------------------------------- */

  ghost->prim.rho = real->prim.rho;
  ghost->prim.u[0] = real->prim.u[0];
  ghost->prim.u[1] = real->prim.u[1];
  ghost->prim.p = real->prim.p;
}






void cell_copy_boundary_data_reflective(cell* real, cell* ghost){
  /* ---------------------------------------------------------
   * Copies the actual data needed for boundaries from a real 
   * cell to a ghost cell. Here for a reflective boundary 
   * condition, where we need to invert the velocities.
   * --------------------------------------------------------- */

  printf("BOUNDARY: Copying from %s -> %s\n", cell_get_index_string(real), cell_get_index_string(ghost));

  ghost->prim.rho = real->prim.rho;
  ghost->prim.u[0] = - real->prim.u[0];
  ghost->prim.u[1] = - real->prim.u[1];
  ghost->prim.p = real->prim.p;
}





void cell_reset_fluxes(){
  /* ----------------------------------------------------
   * reset the fluxes to zero in all cells in the grid
   * ---------------------------------------------------- */

  debugmessage("Resetting fluxes to zero.");

#if NDIM == 1
  for (int i = BC; i < pars.nx + BC; i++){
    gas_init_pstate(&(grid[i].flux));
  }
#elif NDIM == 2
  for (int i = BC; i < pars.nx + BC; i++){
    for (int j = BC; j < pars.nx + BC; j++){
      gas_init_pstate(&(grid[i][j].flux));
    }
  }

#endif
}





float cell_get_total_mass(){
  /* -----------------------------------
   * Compute the total "mass" currently
   * on the grid.
   * ----------------------------------- */

  float mtot = 0;

#if NDIM == 1
  for (int i = BC; i < pars.nx + BC; i++){
    mtot += grid[i].prim.rho;
  }
  mtot *= pars.dx;
  return(mtot);
#elif NDIM == 2
  for (int i = BC; i < pars.nx + BC; i++){
    for (int j = BC; j < pars.nx + BC; j++){
      mtot += grid[i][j].prim.rho;
    }
  }
  mtot *= pars.dx * pars.dx;
  return(mtot);
#endif
}





void cell_print_grid(char field[4]){
  /* -----------------------------------------------
   * Print out some quantity of the entire grid
   * see cell_print_grid_part for options on field[]
   * ----------------------------------------------- */

  int limits[4] = {0, pars.nxtot, 0, pars.nxtot};
  cell_print_grid_part(field, limits);
}





void cell_print_grid_part(char field[4], int* limits){
  /* ---------------------------------------------------------
   * Print out the grid quantities. Select which
   * quantity by field:
   *  "ids": cell ID
   *  "pos": cell position
   *  "rho": density
   *  "v_x": fluid velocity u_x
   *  "v_y": fluid velocity u_x
   *  "vsq": square of fluid velocity u_x^2 + uy^2
   *  "pre": pressure
   *  "frh": density flux
   *  "fux": velocity / momentum flux in x
   *  "fuy": velocity / momentum flux in y
   *  "fpr": pressure flux
   *
   *  limits: array of indices for the boundary of
   *          what to print; i.e. [imin, imax, jmin, jmax]
   *          values at grid indices imin, jmin are 
   *          included in the printout; imax, jmax are not.
   * -------------------------------------------------------- */



#if NDIM == 1
  int imin = limits[0];
  int imax = limits[1];
  for (int i = imin; i < imax; i++){

    if (i==BC || i == pars.nx+BC) printf("|");

    if (strcmp(field, "pos") == 0){
      printf("%8.3f", grid[i].x);
    } else if (strcmp(field, "ids") == 0){
      printf("%8s ", cell_get_index_string(&grid[i]));
    } else if (strcmp(field, "rho") == 0){
      printf("%8.3f", grid[i].prim.rho);
    } else if (strcmp(field, "v_x") == 0){
      printf("%8.3f", grid[i].prim.u[0]);
    } else if (strcmp(field, "v_y") == 0){
      printf("%8.3f", grid[i].prim.u[1]);
    } else if (strcmp(field, "vsq") == 0){
      printf("%8.3f", grid[i].prim.u[0]*grid[i].prim.u[0] + grid[i].prim.u[1]*grid[i].prim.u[1]);
    } else if (strcmp(field, "pre") == 0){
      printf("%8.3f", grid[i].prim.p);
    } else if (strcmp(field, "frh") == 0){
      printf("%8.3f", grid[i].flux.rho);
    } else if (strcmp(field, "fux") == 0){
      printf("%8.3f", grid[i].flux.u[0]);
    } else if (strcmp(field, "fuy") == 0){
      printf("%8.3f", grid[i].flux.u[1]);
    } else if (strcmp(field, "fpr") == 0){
      printf("%8.3f", grid[i].flux.p);
    }
  }
  printf("\n");

#elif NDIM == 2

  int imin = limits[0];
  int imax = limits[1];
  int jmin = limits[2];
  int jmax = limits[3];
  for (int j = jmax - 1; j >= jmin; j--){ /* print highest y first */

    /* print boundary */
    if (j==BC-1 || j == pars.nx+BC-1){
      
      int len = 0;
      if (strcmp(field, "pos") == 0){
        len = 17;
      } else if (strcmp(field, "ids") == 0){
        len = 12;
      } else if (strcmp(field, "rho") == 0){
        len = 8;
      } else if (strcmp(field, "v_x") == 0){
        len = 8;
      } else if (strcmp(field, "v_y") == 0){
        len = 8;
      } else if (strcmp(field, "vsq") == 0){
        len = 8;
      } else if (strcmp(field, "pre") == 0){
        len = 8;
      } else if (strcmp(field, "frh") == 0){
        len = 8;
      } else if (strcmp(field, "fux") == 0){
        len = 8;
      } else if (strcmp(field, "fuy") == 0){
        len = 8;
      } else if (strcmp(field, "fpr") == 0){
        len = 8;
      }

      int dashes = (imax-imin) * len;
      if (imin <= BC) dashes += 2; /* for x boundary */
      if (imax >= pars.nx + BC) dashes += 2; /* for x boundary */
      /* int dashes = pars.nxtot * len + 4; [> +4 for x boundary <] */

      for (int k = 0; k<dashes; k++){ printf("-"); }
      printf("\n");
    }


    for (int i = imin; i < imax; i++){

      /* print boundary */
      if (i==BC || i == pars.nx+BC) printf(" |");

      if (strcmp(field, "pos") == 0){
        printf("(%6.3f, %6.3f) ", grid[i][j].x, grid[i][j].y);
      } else if (strcmp(field, "ids") == 0){
        printf("%12s", cell_get_index_string(&grid[i][j]));
      } else if (strcmp(field, "rho") == 0){
        printf("%8.3f", grid[i][j].prim.rho);
      } else if (strcmp(field, "v_x") == 0){
        printf("%8.3f", grid[i][j].prim.u[0]);
      } else if (strcmp(field, "v_y") == 0){
        printf("%8.3f", grid[i][j].prim.u[1]);
      } else if (strcmp(field, "vsq") == 0){
        printf("%8.3f", grid[i][j].prim.u[0]*grid[i][j].prim.u[0] + grid[i][j].prim.u[1]*grid[i][j].prim.u[1]);
      } else if (strcmp(field, "pre") == 0){
        printf("%8.3f", grid[i][j].prim.p);
      } else if (strcmp(field, "frh") == 0){
        printf("%8.3f", grid[i][j].flux.rho);
      } else if (strcmp(field, "fu[0]") == 0){
        printf("%8.3f", grid[i][j].flux.u[0]);
      } else if (strcmp(field, "fuy") == 0){
        printf("%8.3f", grid[i][j].flux.u[1]);
      } else if (strcmp(field, "fpr") == 0){
        printf("%8.3f", grid[i][j].flux.p);
      }
    }
    printf("\n");
  }

#endif
}




void cell_get_ij(cell* c, int* i, int* j){
  /* ------------------------------------------------ */
  /* Compute the i and j value of a cell such that */
  /* it can be addressed in the grid[] array */
  /* ------------------------------------------------ */

#if NDIM == 1

  *i = c->id;
  *j = 0;

#elif NDIM == 2

  *j = c->id / (pars.nx+2*BC);
  *i = c->id - *j * (pars.nx+2*BC);

#endif
}





char* cell_get_index_string(cell* c){
  /* --------------------------------------------------------
   * get a string for the cell index. Gives a tuple if code 
   * is in 2D. Also correct for ghost boundary cells.
   * if left/lower boundary: returns negative number in a
   * string. If in upper/right boundary: return (nx+int)
   * ------------------------------------------------------- */

  char *output = malloc(MAX_LINE_SIZE * sizeof(char));

#if NDIM == 1
  int id = c->id - BC;
  if (id >= pars.nx){
    sprintf(output, "nx+%d", id - pars.nx);
  } else {
    sprintf(output, "%d", id);
  }

#elif NDIM == 2
  
  strcpy(output, "(");
  char temp[20] = "";
  int i, j;

  j = c->id / pars.nxtot;
  i = c->id - j * pars.nxtot - BC;
  j -= BC; /* only subtract now, otherwise i will be utter bs */

  if (i >= pars.nx){
    sprintf(temp, "nx+%d,", i - pars.nx);
  } else {
    sprintf(temp, "%d,", i);
  }
  strcat(output, temp);

  if (j >= pars.nx){
    sprintf(temp, "nx+%d)", j - pars.nx);
  } else {
    sprintf(temp, "%d)", j);
  }
  strcat(output, temp);

#endif

  return(output);
}
