/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_drag.h - projectile with quadratic air drag: v' = -g y^ - k |v| v.
 * State (x, y, vx, vy); use with pil_rk4_step (n = 4) */
#ifndef PIL_DRAG_H
#define PIL_DRAG_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double g; /* gravity */
  double k; /*  drag coefficient c / m (1/length) */

} pil_drag;

void pil_drag_deriv(double t, const double *y, double *dydt, void *user);

/* Terminal speed sqrt(g / k) (infinite for k = 0) */
double pil_drag_terminal(const pil_drag *p);

/* Trace a launch from (0, h0) with speed v0 at angle theta (radians above
 * horizontal) in steps of dt until it comes back to y = 0 (last point is
 * interpolated onto y = 0) or max_pts points are stored
 * pts[i] = (x, y) at  t = i dt
 * Returns number of points; *t_flight (may be NULL) gets flight time */
int pil_drag_trace(const pil_drag *p, double h0, double v0, double theta,
                   double dt, double (*pts)[2], int max_pts, double *t_flight);

/* Horizontal range from y = h0 back to y = 0 (uses pil_drag_trace) */
double pil_drag_range(const pil_drag *p, double h0, double v0, double theta);

#ifdef __cplusplus
}
#endif

#endif
