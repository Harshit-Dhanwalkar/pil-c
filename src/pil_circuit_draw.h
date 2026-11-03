/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_circuit_draw.h - renderer entry point for circuits (internal) */
#ifndef pil_CIRCUIT_DRAW_H
#define pil_CIRCUIT_DRAW_H

#include "pil/pil_draw.h"
#include "pil_circuit_priv.h"

void pil_circuit_draw(cairo_t *cr, const pil_circuit *c, const pil_view *v,
                      const pil_style *st);

#endif
