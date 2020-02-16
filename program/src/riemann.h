/* top level file for the Riemann solver */


/* Written by Mladen Ivkovic, JAN 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef RIEMANN_H
#define RIEMANN_H


#if RIEMANN == EXACT
#include "riemann/riemann-exact.h"
#elif RIEMANN == TRRS
#include "riemann/riemann-trrs.h"
#elif RIEMANN == TSRS
#include "riemann/riemann-tsrs.h"
#elif RIEMANN == HLL
#include "riemann/riemann-hll.h"
#elif RIEMANN == HLLC
#include "riemann/riemann-hllc.h"
#endif


#ifdef RIEMANN_HLL
#define HLL
#endif
#ifdef RIEMANN_HLLC
#define HLL
#endif




#endif
