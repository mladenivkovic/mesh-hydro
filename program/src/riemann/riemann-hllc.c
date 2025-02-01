/* HLLC Riemann Solver */

/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#include "riemann-hllc.h"
#include "gas.h"
#include "riemann.h"
#include "utils.h"

#include <math.h>

void riemann_solve_hllc(pstate *left, pstate *right, cstate *sol, float xovert,
                        int dimension) {
  /* -------------------------------------------------------------------------
   * Solve the Riemann problem posed by a left and right state.
   * The HLLC solver gives you the fluxes directly, not the states.
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * cstate* sol:     cstate where solution (conserved FLUX) will be written
   * float xovert:    x / t, point where solution shall be sampled
   * int dimension:   which fluid velocity dimension to use. 0: x, 1: y
   * -------------------------------------------------------------------------
   */

  if (riemann_has_vacuum(left, right, dimension)) {
    /* the vacuum solver wants a pstate as the argument for the solution.
     * so give him one, and later translate it back to the conserved flux. */
    pstate pstate_sol;
    gas_init_pstate(&pstate_sol);
    riemann_compute_vacuum_solution(left, right, &pstate_sol, xovert,
                                    dimension);
    gas_get_cflux_from_pstate(&pstate_sol, sol, dimension);

  } else {

    float SL = 0;
    float SR = 0;
    float Sstar = 0;
    riemann_compute_wave_speed_estimates(left, right, &SL, &SR, &Sstar,
                                         dimension);
    riemann_sample_hllc_solution(left, right, SL, SR, Sstar, sol, xovert,
                                 dimension);
  }
}

void riemann_solve_hllc_state(pstate *left, pstate *right, pstate *sol,
                              float xovert, int dimension) {
  /* -------------------------------------------------------------------------
   * Solve the Riemann problem posed by a left and right state
   * Return the state at xovert := x/t instead of the flux
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * cstate* sol:     cstate where solution (conserved FLUX) will be written
   * float xovert:    x / t, point where solution shall be sampled
   * int dimension:   which fluid velocity dimension to use. 0: x, 1: y
   * -------------------------------------------------------------------------
   */

  if (riemann_has_vacuum(left, right, dimension)) {
    /* the vacuum solver wants a pstate as the argument for the solution.
     * so give him one, and later translate it back to the conserved flux. */
    riemann_compute_vacuum_solution(left, right, sol, xovert, dimension);
  } else {
    float SL = 0;
    float SR = 0;
    float Sstar = 0;
    riemann_compute_wave_speed_estimates(left, right, &SL, &SR, &Sstar,
                                         dimension);
    riemann_sample_hllc_solution_state(left, right, SL, SR, Sstar, sol, xovert,
                                       dimension);
  }
}

void riemann_compute_wave_speed_estimates(pstate *left, pstate *right,
                                          float *SL, float *SR, float *Sstar,
                                          int dim) {
  /*--------------------------------------------------------------------------------------------------------
   * Get estimates for the left and right HLLC wave speed.
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * float SL:        left wave speed estimate
   * float SR:        right wave speed estimate
   * float Sstar:     contact wave speed estimate
   * int dim:         which fluid velocity dimension to use. 0: x, 1: y
   * -------------------------------------------------------------------------------------------------------
   */

  /* Start by computint the simple primitive variable speed estimate */
  /* --------------------------------------------------------------- */

  float aL = gas_soundspeed(left);
  float aR = gas_soundspeed(right);

  float temp = 0.25 * (left->rho + right->rho) * (aL + aR);
  float PPV =
      0.5 * (left->p + right->p) + 0.5 * (left->u[dim] - right->u[dim]) * temp;

  if (PPV < 0) {
    PPV = SMALLP;
  }

  float pstar = PPV;
  float ustar =
      0.5 * (left->u[dim] + right->u[dim]) + 0.5 * (left->p - right->p) / temp;

#ifdef HLLC_USE_ADAPTIVE_SPEED_ESTIMATE

  /* use the adaptive wave speed estimate             */
  /* ------------------------------------------------ */

  /* find ratio Q = pmax/pmin, where pmax, pmin are pL and pR */
  float pmin = left->p;
  if (pmin > right->p) {
    pmin = right->p;
  }
  float pmax = left->p;
  if (pmax < right->p) {
    pmax = right->p;
  }
  float qmax = pmax / pmin;

  /* if the ratio pmax/pmin isn't too big, and the primitive variable pressure
   * is between left and right pressure, then PPV approximation is fine */
  if (qmax <= 2. && (pmin <= PPV && PPV <= pmax)) {
    pstar = PPV;
    ustar = 0.5 * (left->u[dim] + right->u[dim]) +
            0.5 * (left->p - right->p) / temp;
  } else {

    if (PPV <= pmin) {
      /* Primitive variable approximation isn't good enough. */
      /* if we expect rarefaction, use the TRRS solver */

      float aLinv = 1. / aL;
      float aRinv = 1. / aR;
      float pLRbeta = powf(left->p / right->p, BETA);

      ustar = ((pLRbeta - 1.) / GM1HALF + left->u[dim] * aLinv * pLRbeta +
               right->u[dim] * aRinv) /
              (aRinv + aLinv * pLRbeta);

      pstar = 0.5 *
              (right->p * pow((1. + aRinv * GM1HALF * (ustar - right->u[dim])),
                              1. / BETA) +
               left->p * pow((1. + aLinv * GM1HALF * (left->u[dim] - ustar)),
                             1. / BETA));
    }

    else {
      /* If not rarefactions, you'll encounter shocks, so use TSRS solver */

      float AL = 2. / GP1 / left->rho;
      float AR = 2. / GP1 / right->rho;
      float BL = GM1OGP1 * left->p;
      float BR = GM1OGP1 * right->p;

      float gL = sqrtf(AL / (PPV + BL));
      float gR = sqrtf(AR / (PPV + BR));

      pstar = (gL * left->p + gR * right->p - (right->u[dim] - left->u[dim])) /
              (gL + gR);
      ustar = 0.5 * (right->u[dim] + left->u[dim] + (pstar - right->p) * gR -
                     (pstar - left->p) * gL);
    }
  }
#endif /* adaptive solution */

  *SL = left->u[dim] - aL * qLR(pstar, left->p);
  *SR = right->u[dim] + aR * qLR(pstar, right->p);
  *Sstar = ustar;
}

float qLR(float pstar, float pLR) {
  /*-----------------------------------------------------
   * Compute q_{L,R} needed for the wave speed estimate.
   *
   * pstar:   (estimated) pressure of the star state
   * pLR:     left or right pressure, depending whether
   *          you want q_L or q_R
   * ------------------------------------------------- */

  if (pstar > pLR) {
    /* shock relation */
    return (sqrtf(1. + 0.5 * (GAMMA + 1.) / GAMMA * (pstar / pLR - 1.)));
  }
  /* Else: rarefaction relation */
  return (1.);
}

void riemann_hllc_compute_star_cstates(pstate *left, pstate *right, float SL,
                                       float SR, float Sstar, cstate *UstarL,
                                       cstate *UstarR, int dim) {

  /*----------------------------------------------------------------------------
   * Compute the !conserved! star states of the solution of the riemann problem
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * float SL:        left wave speed estimate
   * float SR:        right wave speed estimate
   * float Sstar:     contact wave speed estimate
   * cstate* UstarL:  left star conserved state (will be written to)
   * cstate* UstarR:  right star conserved state (will be written to)
   * int dim:         along which dimension are we working
   *------------------------------------------------------------------------- */

  float SLMUL = SL - left->u[dim];
  float SRMUR = SR - right->u[dim];

  /* -------------------------------------------- */
  /* compute left and right conserved star states */
  /* -------------------------------------------- */

  float lcomp = left->rho * SLMUL / (SL - Sstar);
  UstarL->rho = lcomp;
  UstarL->rhou[dim] = lcomp * Sstar;
  UstarL->rhou[(dim + 1) % 2] = lcomp * left->u[(dim + 1) % 2];

  float EL =
      0.5 * left->rho * (left->u[0] * left->u[0] + left->u[1] * left->u[1]) +
      left->p / GM1;
  UstarL->E =
      lcomp * ((EL / left->rho) + (Sstar - left->u[dim]) *
                                      (Sstar + left->p / (left->rho * SLMUL)));

  float rcomp = right->rho * SRMUR / (SR - Sstar);
  UstarR->rho = rcomp;
  UstarR->rhou[dim] = rcomp * Sstar;
  UstarR->rhou[(dim + 1) % 2] = rcomp * right->u[(dim + 1) % 2];

  float ER = 0.5 * right->rho *
                 (right->u[0] * right->u[0] + right->u[1] * right->u[1]) +
             right->p / GM1;
  UstarR->E = rcomp * ((ER / right->rho) +
                       (Sstar - right->u[dim]) *
                           (Sstar + right->p / (right->rho * SRMUR)));
}

void riemann_sample_hllc_solution(pstate *left, pstate *right, float SL,
                                  float SR, float Sstar, cstate *sol,
                                  float xovert, int dim) {
  /*--------------------------------------------------------------------------------------------------
   * Compute the solution of the riemann problem at given time t and x,
   *specified as xovert = x/t Directly returns the flux at x/t, not the state
   *like other Riemann solvers!
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * float SL:        left wave speed estimate
   * float SR:        right wave speed estimate
   * float Sstar:     contact wave speed estimate
   * cstate* sol:     cstate where solution (conserved FLUX) will be written
   * float xovert:    x / t, point where solution shall be sampled
   * int dim:         which fluid velocity direction to use. 0: x, 1: y
   *--------------------------------------------------------------------------------------------------*/

  /* -------------------------------------------- */
  /* compute left and right star conserved states */
  /* -------------------------------------------- */
  cstate UL;
  gas_init_cstate(&UL);
  gas_prim_to_cons(left, &UL);

  cstate UR;
  gas_init_cstate(&UR);
  gas_prim_to_cons(right, &UR);

  cstate UstarL;
  gas_init_cstate(&UstarL);
  cstate UstarR;
  gas_init_cstate(&UstarR);
  riemann_hllc_compute_star_cstates(left, right, SL, SR, Sstar, &UstarL,
                                    &UstarR, dim);

  /* ----------------------------------- */
  /* Compute left and right fluxes       */
  /* ----------------------------------- */
  cstate FL;
  gas_init_cstate(&FL);
  gas_get_cflux_from_cstate(&UL, &FL, dim);

  cstate FR;
  gas_init_cstate(&FR);
  gas_get_cflux_from_cstate(&UR, &FR, dim);

  /* ----------------------------------- */
  /* Compute left and right star fluxes  */
  /* ----------------------------------- */

  cstate FstarL;
  gas_init_cstate(&FstarL);

  FstarL.rho = FL.rho + SL * (UstarL.rho - UL.rho);
  FstarL.rhou[0] = FL.rhou[0] + SL * (UstarL.rhou[0] - UL.rhou[0]);
  FstarL.rhou[1] = FL.rhou[1] + SL * (UstarL.rhou[1] - UL.rhou[1]);
  FstarL.E = FL.E + SL * (UstarL.E - UL.E);

  cstate FstarR;
  gas_init_cstate(&FstarR);

  FstarR.rho = FR.rho + SR * (UstarR.rho - UR.rho);
  FstarR.rhou[0] = FR.rhou[0] + SR * (UstarR.rhou[0] - UR.rhou[0]);
  FstarR.rhou[1] = FR.rhou[1] + SR * (UstarR.rhou[1] - UR.rhou[1]);
  FstarR.E = FR.E + SR * (UstarR.E - UR.E);

  /* ---------------------------- */
  /* finally, sample the solution */
  /* ---------------------------- */

  if (xovert <= SL) {
    /* solution is F_L */
    sol->rho = FL.rho;
    sol->rhou[0] = FL.rhou[0];
    sol->rhou[1] = FL.rhou[1];
    sol->E = FL.E;
  } else if (xovert <= Sstar) {
    /* solution is F*_L */
    sol->rho = FstarL.rho;
    sol->rhou[0] = FstarL.rhou[0];
    sol->rhou[1] = FstarL.rhou[1];
    sol->E = FstarL.E;
  } else if (xovert <= SR) {
    /* solution is F*_R */
    sol->rho = FstarR.rho;
    sol->rhou[0] = FstarR.rhou[0];
    sol->rhou[1] = FstarR.rhou[1];
    sol->E = FstarR.E;
  } else {
    /* solution is F_R */
    sol->rho = FR.rho;
    sol->rhou[0] = FR.rhou[0];
    sol->rhou[1] = FR.rhou[1];
    sol->E = FR.E;
  }
}

void riemann_sample_hllc_solution_state(pstate *left, pstate *right, float SL,
                                        float SR, float Sstar, pstate *sol,
                                        float xovert, int dim) {
  /*--------------------------------------------------------------------------------------------------
   * Compute the solution of the riemann problem at given time t and x,
   *specified as xovert = x/t Returns the state at x/t, not the flux
   *
   * pstate* left:    left state of Riemann problem
   * pstate* right:   right state of Riemann problem
   * float SL:        left wave speed estimate
   * float SR:        right wave speed estimate
   * float Sstar:     contact wave speed estimate
   * cstate* sol:     cstate where solution (conserved FLUX) will be written
   * float xovert:    x / t, point where solution shall be sampled
   * int dim:         which fluid velocity direction to use. 0: x, 1: y
   *--------------------------------------------------------------------------------------------------*/

  /* compute left and right star conserved states */
  cstate UL;
  gas_init_cstate(&UL);
  gas_prim_to_cons(left, &UL);

  cstate UR;
  gas_init_cstate(&UR);
  gas_prim_to_cons(right, &UR);

  cstate UstarL;
  gas_init_cstate(&UstarL);
  cstate UstarR;
  gas_init_cstate(&UstarR);
  riemann_hllc_compute_star_cstates(left, right, SL, SR, Sstar, &UstarL,
                                    &UstarR, dim);

  /* Get left and right primitive star states */
  pstate WstarL;
  gas_init_pstate(&WstarL);
  gas_cons_to_prim(&UstarL, &WstarL);

  pstate WstarR;
  gas_init_pstate(&WstarR);
  gas_cons_to_prim(&UstarR, &WstarR);

  /* finally, sample the solution */
  if (xovert <= SL) {
    /* solution is W_L */
    sol->rho = left->rho;
    sol->u[0] = left->u[0];
    sol->u[1] = left->u[1];
    sol->p = left->p;
  } else if (xovert <= Sstar) {
    /* solution is W*_L */
    sol->rho = WstarL.rho;
    sol->u[0] = WstarL.u[0];
    sol->u[1] = WstarL.u[1];
    sol->p = WstarL.p;
  } else if (xovert <= SR) {
    /* solution is W*_R */
    sol->rho = WstarR.rho;
    sol->u[0] = WstarR.u[0];
    sol->u[1] = WstarR.u[1];
    sol->p = WstarR.p;
  } else {
    /* solution is W_R */
    sol->rho = right->rho;
    sol->u[0] = right->u[0];
    sol->u[1] = right->u[1];
    sol->p = right->p;
  }
}

void riemann_get_hllc_full_solution_for_WAF(pstate *left, pstate *right,
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
   * F*_R, F_R float delta_q[3]:  differences in densities over all four waves:
   *                      U*_L - U_L, U*_R - U*_L, U_R - U*_R
   * int dim:           which fluid velocity direction to use. 0: x, 1: y
   * ------------------------------------------------------------------------------------------
   */

  /* first compute wave speeds */

  float SL = 0;
  float SR = 0;
  float Sstar = 0;
  riemann_compute_wave_speed_estimates(left, right, &SL, &SR, &Sstar, dim);

  /* store wave speeds */
  S[0] = SL;
  S[1] = Sstar;
  S[2] = SR;

  /* --------------------------------------- */
  /* compute left and right conserved states */
  /* --------------------------------------- */
  cstate UL;
  gas_init_cstate(&UL);
  gas_prim_to_cons(left, &UL);

  cstate UR;
  gas_init_cstate(&UR);
  gas_prim_to_cons(right, &UR);

  /* -------------------------------------------- */
  /* compute left and right conserved star states */
  /* -------------------------------------------- */

  /* compute left and right star conserved states */
  cstate UstarL;
  gas_init_cstate(&UstarL);
  cstate UstarR;
  gas_init_cstate(&UstarR);
  riemann_hllc_compute_star_cstates(left, right, SL, SR, Sstar, &UstarL,
                                    &UstarR, dim);

  /* store jumps over density in delta_q */
  delta_q[0] = UstarL.rho - UL.rho;
  delta_q[1] = UstarR.rho - UstarL.rho;
  delta_q[2] = UR.rho - UstarR.rho;

  /* ----------------------------------- */
  /* Compute left and right fluxes       */
  /* ----------------------------------- */
  cstate FL;
  gas_init_cstate(&FL);
  gas_get_cflux_from_cstate(&UL, &FL, dim);

  cstate FR;
  gas_init_cstate(&FR);
  gas_get_cflux_from_cstate(&UR, &FR, dim);

  /* ----------------------------------- */
  /* Compute left and right star fluxes  */
  /* ----------------------------------- */

  cstate FstarL;
  gas_init_cstate(&FstarL);

  FstarL.rho = FL.rho + SL * (UstarL.rho - UL.rho);
  FstarL.rhou[0] = FL.rhou[0] + SL * (UstarL.rhou[0] - UL.rhou[0]);
  FstarL.rhou[1] = FL.rhou[1] + SL * (UstarL.rhou[1] - UL.rhou[1]);
  FstarL.E = FL.E + SL * (UstarL.E - UL.E);

  cstate FstarR;
  gas_init_cstate(&FstarR);

  FstarR.rho = FR.rho + SR * (UstarR.rho - UR.rho);
  FstarR.rhou[0] = FR.rhou[0] + SR * (UstarR.rhou[0] - UR.rhou[0]);
  FstarR.rhou[1] = FR.rhou[1] + SR * (UstarR.rhou[1] - UR.rhou[1]);
  FstarR.E = FR.E + SR * (UstarR.E - UR.E);

  /* store the fluxes */
  fluxes[0] = FL;
  fluxes[1] = FstarL;
  fluxes[2] = FstarR;
  fluxes[3] = FR;
}

void riemann_compute_star_states(pstate *left, pstate *right, float *pstar,
                                 float *ustar, int dimension) {
  /* The WAF scheme requires this function when calling
   * riemann_get_full_solution, but it's done differently
   * for the HLLC solver (riemann_get_hllc_full_solution)
   * so this function is just empty here so that the code
   * will compile correctly. */
  throw_error(
      "riemann_compute_star_states shouldn't be called for the HLLC solver");
}
