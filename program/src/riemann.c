/* functions that all Riemann solvers use */

/* Written by Mladen Ivkovic, FEB 2020
 * mladen.ivkovic@hotmail.com           */

#include <math.h>
#include <stddef.h>

#include "defines.h"
#include "gas.h"
#include "params.h"
#include "riemann.h"

#if RIEMANN == HLLC
#include "utils.h"
#endif

extern params pars;

void riemann_solve(pstate *left, pstate *right, pstate *sol, float xovert,
                   int dimension) {
  /* -------------------------------------------------------------------------
   * Solve the Riemann problem posed by a left and right state
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * pstate* sol:     pstate where solution will be written
   * float xovert:    x / t, point where solution shall be sampled
   * int dimension:   which fluid velocity dimension to use. 0: x, 1: y
   * -------------------------------------------------------------------------
   */

  if (riemann_has_vacuum(left, right, dimension)) {
    riemann_compute_vacuum_solution(left, right, sol, xovert, dimension);
  } else {
    float pstar = 0;
    float ustar = 0;
    riemann_compute_star_states(left, right, &pstar, &ustar, dimension);
    riemann_sample_solution(left, right, pstar, ustar, sol, xovert, dimension);
  }
}

int riemann_has_vacuum(pstate *left, pstate *right, int dimension) {
  /* -------------------------------------------------------------------------
   * Check whether we work with vacuum
   * returns true (1) if vacuum, 0 otherwise
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * int dimension:   which fluid velocity dimension to use. 0: x, 1: y
   * -------------------------------------------------------------------------
   */

  if (left->rho <= SMALLRHO) {
    return (1);
  }
  if (right->rho <= SMALLRHO) {
    return (1);
  }

  float delta_u = right->u[dimension] - left->u[dimension];
  float u_crit = 2. / GM1 * (gas_soundspeed(left) + gas_soundspeed(right));

  if (delta_u < u_crit) {
    return (0);
  }
  return (1);

}

void riemann_compute_vacuum_solution(pstate *left, pstate *right, pstate *sol,
                                     float xovert, int dim) {
  /* -------------------------------------------------------------------------
   * Solve the Riemann problem posed by a left and right state
   * and sample the solution at given x and t
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * pstate* sol:     pstate where solution will be written
   * float xovert:    x / t, point where solution shall be sampled
   * int dim:   which fluid velocity dim to use. 0: x, 1: y
   * -------------------------------------------------------------------------
   */

  int notdim = (dim + 1) % 2;

  if (left->rho <= SMALLRHO && right->rho <= SMALLRHO) {
    sol->rho = SMALLRHO;
    sol->u[dim] = SMALLU;
    sol->u[notdim] = SMALLU;
    sol->p = SMALLP;
    return;
  }

  if (left->rho <= SMALLRHO) {
    /*------------------------*/
    /* Left vacuum state      */
    /*------------------------*/
    float aR = gas_soundspeed(right);
    float SR = right->u[dim] - 2. * aR / GM1; /* vacuum front speed */
    float SHR = right->u[dim] + aR; /* speed of head of right rarefaction fan */

    if (xovert <= SR) {
      /* left vacuum */
      sol->rho = SMALLRHO;
#ifdef USE_AS_RIEMANN_SOLVER
      sol->u[dim] = SR;
#else
      sol->u[dim] = SMALLU;
#endif
      sol->u[notdim] = SMALLU;
      sol->p = SMALLP;
    } else if (xovert < SHR) {
      /* inside rarefaction */
      float precomp =
          pow((2. / GP1 - GM1OGP1 / aR * (right->u[dim] - xovert)), (2. / GM1));
      sol->rho = right->rho * precomp;
      sol->u[dim] = 2. / GP1 * (GM1HALF * right->u[dim] - aR + xovert);
      sol->u[notdim] = right->u[notdim];
      sol->p = right->p * pow(precomp, GAMMA);
    } else {
      /* original right pstate */
      sol->rho = right->rho;
      sol->u[dim] = right->u[dim];
      sol->u[notdim] = right->u[notdim];
      sol->p = right->p;
    }
  }

  else if (right->rho <= SMALLRHO) {
    /*------------------------*/
    /* Right vacuum state     */
    /*------------------------*/

    float aL = gas_soundspeed(left);
    float SL = left->u[dim] + 2. * aL / GM1; /* vacuum front speed */
    float SHL = left->u[dim] - aL; /* speed of head of left rarefaction fan */

    if (xovert >= SL) {
      /* right vacuum */
      sol->rho = SMALLRHO;
#ifdef USE_AS_RIEMANN_SOLVER
      sol->u[dim] = SL;
#else
      sol->u[dim] = SMALLU;
#endif
      sol->u[notdim] = SMALLU;
      sol->p = SMALLP;
    } else if (xovert > SHL) {
      /* inside rarefaction */
      float precomp =
          pow((2. / GP1 + GM1OGP1 / aL * (left->u[dim] - xovert)), (2. / GM1));
      sol->rho = left->rho * precomp;
      sol->u[dim] = 2. / GP1 * (GM1HALF * left->u[dim] + aL + xovert);
      sol->u[notdim] = left->u[notdim];
      sol->p = left->p * pow(precomp, GAMMA);
    } else {
      /* original left pstate */
      sol->rho = left->rho;
      sol->u[dim] = left->u[dim];
      sol->u[notdim] = left->u[notdim];
      sol->p = left->p;
    }
  } else {
    /*------------------------*/
    /* Vacuum generating case */
    /*------------------------*/

    float aL = gas_soundspeed(left);
    float aR = gas_soundspeed(right);
    float SL = left->u[dim] + 2. * aL / GM1;  /* vacuum front speed */
    float SR = right->u[dim] - 2. * aR / GM1; /* vacuum front speed */
    float SHL = left->u[dim] - aL;  /* speed of head of left rarefaction fan */
    float SHR = right->u[dim] + aR; /* speed of head of right rarefaction fan */

    if (xovert <= SHL) {
      /* left original pstate*/
      sol->rho = left->rho;
      sol->u[dim] = left->u[dim];
      sol->u[notdim] = left->u[notdim];
      sol->p = left->p;
    } else if (xovert < SL) {
      /* inside rarefaction fan from right to left */
      float precomp =
          pow((2. / GP1 + GM1OGP1 / aL * (left->u[dim] - xovert)), (2. / GM1));
      sol->rho = left->rho * precomp;
      sol->u[dim] = 2. / GP1 * (GM1HALF * left->u[dim] + aL + xovert);
      sol->u[notdim] = left->u[notdim];
      sol->p = left->p * pow(precomp, GAMMA);
    } else if (xovert < SR) {
      /* vacuum region */
      sol->rho = SMALLRHO;
#ifdef USE_AS_RIEMANN_SOLVER
      sol->u[dim] = 0.5 * (SL + SR);
#else
      sol->u[dim] = SMALLU;
#endif
      sol->u[notdim] = SMALLU;
      sol->p = SMALLP;
    } else if (xovert < SHR) {
      /* inside rarefaction fan from left to right */
      float precomp =
          pow((2. / GP1 - GM1OGP1 / aR * (right->u[dim] - xovert)), (2. / GM1));
      sol->rho = right->rho * precomp;
      sol->u[dim] = 2. / GP1 * (GM1HALF * right->u[dim] - aR + xovert);
      sol->u[notdim] = right->u[notdim];
      sol->p = right->p * pow(precomp, GAMMA);
    } else {
      /* right original pstate */
      sol->rho = right->rho;
      sol->u[dim] = right->u[dim];
      sol->u[notdim] = right->u[notdim];
      sol->p = right->p;
    }
  }

}

void riemann_sample_solution(pstate *left, pstate *right, float pstar,
                             float ustar, pstate *sol, float xovert, int dim) {
  /*--------------------------------------------------------------------------------------------------
   * Compute the solution of the riemann problem at given time t and x,
   *specified as xovert = x/t
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * float pstar:     pressure of star region
   * float ustar:     velocity of star region
   * pstate* sol:     pstate where solution will be written
   * float xovert:    x / t, point where solution shall be sampled
   * int dim:         which fluid velocity direction to use. 0: x, 1: y
   *--------------------------------------------------------------------------------------------------*/

  int otherdim = (dim + 1) % 2;

  if (xovert <= ustar) {
    /*------------------------*/
    /* We're on the left side */
    /*------------------------*/
    float aL = gas_soundspeed(left);
    float pstaroverpL = pstar / left->p;

    if (pstar <= left->p) {
      /*------------------*/
      /* left rarefaction */
      /*------------------*/
      float SHL = left->u[dim] - aL; /* speed of head of left rarefaction fan */
      if (xovert < SHL) {
        /* we're outside the rarefaction fan */
        sol->rho = left->rho;
        sol->u[dim] = left->u[dim];
        sol->u[otherdim] = left->u[otherdim];
        sol->p = left->p;
      } else {
        float astarL = aL * powf(pstaroverpL, BETA);
        float STL = ustar - astarL; /* speed of tail of left rarefaction fan */
        if (xovert < STL) {
          /* we're inside the fan */
          float precomp = pow(
              (2. / GP1 + GM1OGP1 / aL * (left->u[dim] - xovert)), (2. / GM1));
          sol->rho = left->rho * precomp;
          sol->u[dim] = 2. / GP1 * (GM1HALF * left->u[dim] + aL + xovert);
          sol->u[otherdim] = left->u[otherdim];
          sol->p = left->p * pow(precomp, GAMMA);
        } else {
          /* we're in the star region */
          sol->rho = left->rho * powf(pstaroverpL, ONEOVERGAMMA);
          sol->u[dim] = ustar;
          sol->u[otherdim] = left->u[otherdim];
          sol->p = pstar;
        }
      }
    } else {
      /*------------------*/
      /* left shock       */
      /*------------------*/
      float SL = left->u[dim] - aL * sqrtf(0.5 * GP1 / GAMMA * pstaroverpL +
                                           BETA); /* left shock speed */
      if (xovert < SL) {
        /* we're outside the shock */
        sol->rho = left->rho;
        sol->u[dim] = left->u[dim];
        sol->u[otherdim] = left->u[otherdim];
        sol->p = left->p;
      } else {
        /* we're in the star region */
        sol->rho =
            (pstaroverpL + GM1OGP1) / (GM1OGP1 * pstaroverpL + 1.) * left->rho;
        sol->u[dim] = ustar;
        sol->u[otherdim] = left->u[otherdim];
        sol->p = pstar;
      }
    }
  } else {
    /*-------------------------*/
    /* We're on the right side */
    /*-------------------------*/
    float aR = gas_soundspeed(right);
    float pstaroverpR = pstar / right->p;
    if (pstar <= right->p) {

      /*-------------------*/
      /* right rarefaction */
      /*-------------------*/
      float SHR =
          right->u[dim] + aR; /* speed of head of right rarefaction fan */
      if (xovert > SHR) {
        /* we're outside the rarefaction fan */
        sol->rho = right->rho;
        sol->u[dim] = right->u[dim];
        sol->u[otherdim] = right->u[otherdim];
        sol->p = right->p;
      } else {
        float astarR = aR * powf(pstaroverpR, BETA);
        float STR = ustar + astarR; /* speed of tail of right rarefaction fan */
        if (xovert > STR) {
          /* we're inside the fan */
          float precomp = pow(
              (2. / GP1 - GM1OGP1 / aR * (right->u[dim] - xovert)), (2. / GM1));
          sol->rho = right->rho * precomp;
          sol->u[dim] = 2. / GP1 * (GM1HALF * right->u[dim] - aR + xovert);
          sol->u[otherdim] = right->u[otherdim];
          sol->p = right->p * pow(precomp, GAMMA);
        } else {
          /* we're in the star region */
          sol->rho = right->rho * powf(pstaroverpR, ONEOVERGAMMA);
          sol->u[dim] = ustar;
          sol->u[otherdim] = right->u[otherdim];
          sol->p = pstar;
        }
      }
    } else {
      /*------------------*/
      /* right shock      */
      /*------------------*/
      float SR = right->u[dim] + aR * sqrtf(0.5 * GP1 / GAMMA * pstaroverpR +
                                            BETA); /* right shock speed */
      if (xovert > SR) {
        /* we're outside the shock */
        sol->rho = right->rho;
        sol->u[dim] = right->u[dim];
        sol->u[otherdim] = right->u[otherdim];
        sol->p = right->p;
      } else {
        /* we're in the star region */
        sol->rho =
            (pstaroverpR + GM1OGP1) / (GM1OGP1 * pstaroverpR + 1.) * right->rho;
        sol->u[dim] = ustar;
        sol->u[otherdim] = right->u[otherdim];
        sol->p = pstar;
      }
    }
  }

}

void riemann_get_full_solution_for_WAF(pstate *left, pstate *right, float S[3],
                                       cstate fluxes[4], float delta_q[3],
                                       int dim) {
  /*-------------------------------------------------------------------------------------------
   * Compute (and "return") the full solution of the Riemann problem: Get all
   * wave speeds, the fluxes of all four states U_L, U*_L, U*_R, U_R, and the
   * difference in densities between each wave. This function is needed for the
   * WAF method, where we sum up all the occuring fluxes with different weights.
   * This function handles the vacuum case.
   *
   * pstate* left:      left primitive state of Riemann problem
   * pstate* right:     right primitive state of Riemann problem
   * float S[3]:        where wave speeds will be written to
   * cstate fluxes[4]:  where the four fluxes will be written to: F_L, F*_L,
   * F*_R, F_R float delta_q[3]:  differences in densities over all four waves:
   *                      U*_L - U_L, U*_R - U*_L, U_R - U*_R
   * int dim:           which fluid velocity direction to use. 0: x, 1: y
   * ------------------------------------------------------------------------------------------
   */

#if RIEMANN == HLLC
  throw_error("In riemann_get_full_solution: You shouldn't be calling this "
              "function with the HLLC Riemann solver.");
#endif

  int notdim = (dim + 1) % 2;

  /* First, get the star velocity and pressure */
  float pstar, ustar;
  riemann_compute_star_states(left, right, &pstar, &ustar, dim);

  /* Now fill them into left and right star state variables */
  pstate star_left;
  gas_init_pstate(&star_left);
  star_left.u[dim] = ustar;
  star_left.u[notdim] = left->u[notdim];
  star_left.p = pstar;
  float pstaroverpL = pstar / left->p;
  if (pstar <= left->p) {
    /* rarefaction star state */
    star_left.rho = left->rho * powf(pstaroverpL, ONEOVERGAMMA);
  } else {
    /* shock star state */
    star_left.rho =
        (pstaroverpL + GM1OGP1) / (GM1OGP1 * pstaroverpL + 1.) * left->rho;
  }

  pstate star_right;
  gas_init_pstate(&star_right);
  star_right.u[dim] = ustar;
  star_right.u[notdim] = right->u[notdim];
  star_right.p = pstar;
  float pstaroverpR = pstar / right->p;
  if (pstar <= right->p) {
    /* rarefaction star state */
    star_right.rho = right->rho * powf(pstaroverpR, ONEOVERGAMMA);
  } else {
    /* shock star state */
    star_right.rho =
        (pstaroverpR + GM1OGP1) / (GM1OGP1 * pstaroverpR + 1.) * right->rho;
  }

  /* get wave speeds. We always expect 3 waves. */

  float aL = gas_soundspeed(left);
  if (pstar <= left->p) {
    /* left rarefaction head velocity */
    S[0] = (left->u[dim] - aL);
  } else {
    /* left shock velocity */
    S[0] = (left->u[dim] - aL * sqrtf(0.5 * GP1 / GAMMA * pstaroverpL + BETA));
  }

  S[1] = ustar;

  float aR = gas_soundspeed(right);
  if (pstar <= right->p) {
    /* right rarefaction head velocity */
    S[2] = (right->u[dim] + aR);
  } else {
    /* right shock velocity */
    S[2] = (right->u[dim] + aR * sqrtf(0.5 * GP1 / GAMMA * pstaroverpR + BETA));
  }

  /* find solution at x = 0 in case we need it */
  pstate state_at_x_zero;
  gas_init_pstate(&state_at_x_zero);
  riemann_sample_solution(left, right, pstar, ustar, &state_at_x_zero,
                          /*xovert =*/0., dim);

  /* get fluxes from the states. For 1D Euler equations, there will always be 4
   * fluxes to consider. */

  /* state 1 */
  gas_get_cflux_from_pstate(left, &fluxes[0], dim);

  /* state 2 */
  if (pstar <= left->p) {
    /* rarefaction. But is it sonic? */
    /* use approximation for fan instead of adding one more wave and properly
     * integrating */
    float astarL = gas_soundspeed(&star_left);
    float tailL = ustar - astarL;
    if (tailL * S[0] < 0.) {
      /* if head and tail speeds have different sign, we have a sonic
       * rarefaction. */
      gas_get_cflux_from_pstate(&state_at_x_zero, &fluxes[1], dim);
    } else {
      /* Non-sonic rarefaction */
      gas_get_cflux_from_pstate(&star_left, &fluxes[1], dim);
    }
  } else {
    /* Shock */
    gas_get_cflux_from_pstate(&star_left, &fluxes[1], dim);
  }

  /* state 3 */
  if (pstar <= right->p) {
    /* rarefaction. But is it sonic? */
    /* use approximation for fan instead of adding one more wave and properly
     * integrating */
    float astarR = gas_soundspeed(&star_right);
    float tailR = ustar + astarR;
    if (tailR * S[2] < 0.) {
      /* if head and tail speeds have different sign, we have asonic
       * rarefaction. */
      gas_get_cflux_from_pstate(&state_at_x_zero, &fluxes[2], dim);
    } else {
      /* Non-sonic rarefaction */
      gas_get_cflux_from_pstate(&star_right, &fluxes[2], dim);
    }
  } else {
    /* shock */
    gas_get_cflux_from_pstate(&star_right, &fluxes[2], dim);
  }

  /* state 4 */
  gas_get_cflux_from_pstate(right, &fluxes[3], dim);

  delta_q[0] = star_left.rho - left->rho;
  delta_q[1] = star_right.rho - star_left.rho;
  delta_q[2] = right->rho - star_right.rho;
}

void riemann_get_full_vacuum_solution_for_WAF(pstate *left, pstate *right,
                                              float S[3], cstate fluxes[4],
                                              float delta_q[3], int dim) {
  /*-------------------------------------------------------------------------------------------
   * Compute (and "return") the full solution of the Riemann problem: Get all
   * wave speeds, the fluxes of all four states U_L, U*_L, U*_R, U_R, and the
   * difference in densities between each wave. This function is needed for the
   * WAF method, where we sum up all the occuring fluxes with different weights.
   * This function handles the vacuum case.
   *
   * pstate* left:      left primitive state of Riemann problem
   * pstate* right:     right primitive state of Riemann problem
   * float S[3]:        where wave speeds will be written to
   * cstate fluxes[4]:  where the four fluxes will be written to: F_L, F*_L,
   *                    F*_R, F_R
   * float delta_q[3]:  differences in densities over all four waves:
   *                      U*_L - U_L, U*_R - U*_L, U_R - U*_R
   * int dim:           which fluid velocity direction to use. 0: x, 1: y
   * ------------------------------------------------------------------------------------------
   */

  int notdim = (dim + 1) % 2;

  pstate vacuum;
  gas_init_pstate(&vacuum);

  /* Easy case: Both are vacuum */
  if (left->rho <= SMALLRHO && right->rho <= SMALLRHO) {
    for (int k = 0; k < 4; k++) {
      gas_get_cflux_from_pstate(&vacuum, &fluxes[k], dim);
    }
    for (int k = 0; k < 3; k++) {
      S[k] = 0.;
      delta_q[k] = 0.;
    }
    return;
  }

  pstate state_at_x_zero;
  gas_init_pstate(&state_at_x_zero);
  pstate star_left;
  gas_init_pstate(&star_left);
  pstate star_right;
  gas_init_pstate(&star_right);

  if (left->rho <= SMALLRHO) {
    /*------------------------*/
    /* Left vacuum state      */
    /*------------------------*/
    float aR = gas_soundspeed(right);
    float SR = right->u[dim] - 2. * aR / GM1; /* vacuum front speed */
    float SHR = right->u[dim] + aR; /* speed of head of right rarefaction fan */
    float rho_star_right = 0;

    /* We always have a rarefaction, so find the state at x = 0 */
    if (0. <= SR) {
      state_at_x_zero.rho = SMALLRHO;
      state_at_x_zero.u[dim] = SMALLU;
      state_at_x_zero.u[notdim] = SMALLU;
      state_at_x_zero.p = SMALLP;
    } else if (0. < SHR) {
      /* inside rarefaction */
      float precomp = pow(
          (2. / GP1 - GM1OGP1 / aR * (right->u[dim] /*- x/t */)), (2. / GM1));
      state_at_x_zero.rho = right->rho * precomp;
      state_at_x_zero.u[dim] =
          2. / GP1 * (GM1HALF * right->u[dim] - aR /* + x/t */);
      state_at_x_zero.u[notdim] = right->u[notdim];
      state_at_x_zero.p = right->p * pow(precomp, GAMMA);
    } else {
      /* original right pstate */
      state_at_x_zero.rho = right->rho;
      state_at_x_zero.u[dim] = right->u[dim];
      state_at_x_zero.u[notdim] = right->u[notdim];
      state_at_x_zero.p = right->p;
    }

    /* compute the fluxes */
    gas_get_cflux_from_pstate(left, &fluxes[0], dim);
    gas_get_cflux_from_pstate(&vacuum, &fluxes[1], dim);
    /* use approximation for fan instead of adding one more wave and properly
     * integrating */
    if (SR * SHR < 0) {
      /* Sonic rarefaction */
      gas_get_cflux_from_pstate(&state_at_x_zero, &fluxes[2], dim);
      rho_star_right = state_at_x_zero.rho;
    } else {
      /* Non-sonic rarefaction */
      gas_get_cflux_from_pstate(right, &fluxes[2], dim);
      rho_star_right = right->rho;
    }
    gas_get_cflux_from_pstate(right, &fluxes[3], dim);

    /* store the wave speeds */
    S[0] = SMALLU;
    S[1] = SR;
    S[2] = SHR;

    /* get density jumps */
    delta_q[0] = 0.;
    delta_q[1] = 0.;
    delta_q[2] = right->rho - rho_star_right;

  }

  else if (right->rho <= SMALLRHO) {
    /*------------------------*/
    /* Right vacuum state     */
    /*------------------------*/

    float aL = gas_soundspeed(left);
    float SL = left->u[dim] + 2. * aL / GM1; /* vacuum front speed */
    float SHL = left->u[dim] - aL; /* speed of head of left rarefaction fan */
    float rho_star_left = 0;

    if (0. >= SL) {
      /* right vacuum */
      state_at_x_zero.rho = SMALLRHO;
      state_at_x_zero.u[dim] = SMALLU;
      state_at_x_zero.u[notdim] = SMALLU;
      state_at_x_zero.p = SMALLP;
    } else if (0. > SHL) {
      /* inside rarefaction */
      float precomp = pow(
          (2. / GP1 + GM1OGP1 / aL * (left->u[dim] /* - x/t */)), (2. / GM1));
      state_at_x_zero.rho = left->rho * precomp;
      state_at_x_zero.u[dim] =
          2. / GP1 * (GM1HALF * left->u[dim] + aL /* + x/t */);
      state_at_x_zero.u[notdim] = left->u[notdim];
      state_at_x_zero.p = left->p * pow(precomp, GAMMA);
    } else {
      /* original left pstate */
      state_at_x_zero.rho = left->rho;
      state_at_x_zero.u[dim] = left->u[dim];
      state_at_x_zero.u[notdim] = left->u[notdim];
      state_at_x_zero.p = left->p;
    }

    /* compute the fluxes */
    gas_get_cflux_from_pstate(left, &fluxes[0], dim);
    /* use approximation for fan instead of adding one more wave and properly
     * integrating */
    if (SL * SHL < 0) {
      /* Sonic rarefaction */
      gas_get_cflux_from_pstate(&state_at_x_zero, &fluxes[1], dim);
      rho_star_left = state_at_x_zero.rho;
    } else {
      /* Non-sonic rarefaction */
      gas_get_cflux_from_pstate(left, &fluxes[1], dim);
      rho_star_left = left->rho;
    }
    gas_get_cflux_from_pstate(&vacuum, &fluxes[2], dim);
    gas_get_cflux_from_pstate(right, &fluxes[3], dim);

    /* store the wave speeds */
    S[0] = SHL;
    S[1] = SL;
    S[2] = SMALLU;

    /* get density jumps */
    delta_q[0] = rho_star_left - left->rho;
    delta_q[1] = 0.;
    delta_q[2] = 0.;
  }

  else {
    /*------------------------*/
    /* Vacuum generating case */
    /*------------------------*/

    float aL = gas_soundspeed(left);
    float aR = gas_soundspeed(right);
    float SL = left->u[dim] + 2. * aL / GM1;  /* vacuum front speed */
    float SR = right->u[dim] - 2. * aR / GM1; /* vacuum front speed */
    float SHL = left->u[dim] - aL;  /* speed of head of left rarefaction fan */
    float SHR = right->u[dim] + aR; /* speed of head of right rarefaction fan */

    float speed_average, precomp, f;

    pstate star_left;
    gas_init_pstate(&star_left);

    speed_average = 0.5 * (SHL + SL);
    precomp =
        pow((2. / GP1 + GM1OGP1 / aL * (left->u[dim] - /*x/t=*/speed_average)),
            (2. / GM1));
    star_left.rho = left->rho * precomp;
    star_left.u[dim] = 0.5 * (SL + SR);
    star_left.p = left->p * pow(precomp, GAMMA);

    pstate star_right;
    gas_init_pstate(&star_right);
    speed_average = 0.5 * (SHR + SR);
    precomp =
        pow((2. / GP1 - GM1OGP1 / aR * (right->u[dim] - /*x/t=*/speed_average)),
            (2. / GM1));
    star_right.rho = right->rho * precomp;
    star_right.u[dim] = 0.5 * (SL + SR);
    star_right.p = right->p * pow(precomp, GAMMA);

    f = (SHL - SL) / (SHL - 0.5 * (SL + SR));
    star_left.rho *= f;
    star_left.p *= f;

    f = (SHR - SR) / (SHR - 0.5 * (SL + SR));
    star_right.rho *= precomp;
    star_right.p *= precomp;

    /* compute the fluxes */
    gas_get_cflux_from_pstate(left, &fluxes[0], dim);
    gas_get_cflux_from_pstate(&star_left, &fluxes[1], dim);
    gas_get_cflux_from_pstate(&star_right, &fluxes[2], dim);
    gas_get_cflux_from_pstate(right, &fluxes[3], dim);

    /* store the wave speeds */
    S[0] = SHL;
    S[1] = 0.5 * (SL + SR);
    S[2] = SHR;

    /* get density jumps */
    delta_q[0] = star_left.rho - left->rho;
    delta_q[1] = star_right.rho - star_left.rho;
    delta_q[2] = right->rho - star_right.rho;
  }

}
