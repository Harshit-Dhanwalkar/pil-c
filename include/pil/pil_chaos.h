/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_chaos.h - driven damped pendulum and Poincare sections */
#ifndef pil_CHAOS_H
#define pil_CHAOS_H

#ifdef __cplusplus
extern "C" {
#endif

/* theta'' + theta' / q + sin(theta) = g cos(w_d t), state (theta, omega);
 * chaotic for q = 2, w_d = 2/3, g = 1.5. Use with pil_rk4_step (n = 2) */
typedef struct {
  double q, w_d, g;
} pil_driven_pendulum;
void pil_driven_pendulum_deriv(double t, const double *y, double *dydt,
                               void *user);

/* Integrate nper drive periods with steps_per_period RK4 steps each
 * Fills stroboscopic section (theta wrapped to [-pi, pi), omega) at t = k T
 * into section[k] (nper entries) and when path is non-NULL, every step into
 * path[i] (nper * steps_per_period entries) */
void pil_poincare_section(const pil_driven_pendulum *p, double theta0,
                          double omega0, int nper, int steps_per_period,
                          double (*section)[2], double (*path)[2]);

#ifdef __cplusplus
}
#endif

#endif
