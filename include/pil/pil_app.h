/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/*
 * pil_app.h - optional interactive window (SDL2 + Cairo)
 *
 *   pil_app_config cfg = pil_app_defaults("My problem", 900, 480);
 *   pil_view_fit(&cfg.view, -2, -0.2, 2, 1.6, 900, 480, 24);
 *   pil_app_run(world, &cfg);
 *
 * Keys:   Space pause/resume    R reset    Right step (while paused)
 *         + / -  simulation speed          H toggle HUD
 *         P  save PNG   F  save PDF   G  save SVG     Esc quit
 * Mouse:  left drag  pull a body with a soft spring
 *         wheel zoom,  right/middle drag pan
 */
#ifndef pil_APP_H
#define pil_APP_H

#include "pil_core.h"
#include "pil_draw.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  const char *title;
  int width, height;
  pil_view view;
  pil_style style;
  double time_scale; /* 1.0 = real time                                  */
  double max_dt;     /* largest integration step (default 5e-4)          */
  double grab_k;     /* mouse spring stiffness per kg (default 300)      */
  int start_paused;
  int show_hud;
  /* automation / testing: exit after this many frames (0 = never) and save a
   * PNG on exit */
  int max_frames;
  const char *exit_screenshot;
  /* optional per-frame hook, called before drawing with the simulation time */
  void (*on_frame)(pil_world *w, double t, void *user);
  void *user;
} pil_app_config;

pil_app_config pil_app_defaults(const char *title, int width, int height);
int pil_app_run(pil_world *w, pil_app_config *cfg);

#ifdef __cplusplus
}
#endif
#endif
