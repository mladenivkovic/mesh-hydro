/* top level file for sources */

/* Written by Mladen Ivkovic, APR 2020
 * mladen.ivkovic@hotmail.com           */

#ifndef SOURCES_H
#define SOURCES_H

#include "gas.h"

void sources_get_acceleration(void);

void sources_get_source_vector(cstate* s, const float acc[2], cstate* cons);

void sources_update_state(float dt);

#endif
