/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_plot.h - renderer-side entry point for plot panels (internal) */

#ifndef pil_PLOT_H
#define pil_PLOT_H

#include "pil/pil_draw.h"
#include "pil_internal.h"

void pil_plot_draw(cairo_t *cr, const pil_plot *p, const pil_world *w,
                   const pil_view *v, const pil_style *st);

#endif
