/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_kapitza.h - Kapitza pendulum: a pendulum whose pivot oscillates
 * vertically as y_p = a cos(w t)
 * angle theta is measured from upward vertical (theta = 0 is inverted
 * position):
 *   theta'' = (g - a w^2 cos(w t)) / l * sin(theta) - gamma theta'
 * Averaged over fast vibration bob moves in effective potential:
 *   V_eff = -(g / l) cos(theta) + a^2 w^2 / (4 l^2) sin^2(theta)
 * so inverted position is stable when a^2 w^2 > 2 g l */

#ifndef PIL_KAPITZA_H
#define PIL_KAPITZA_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double g, l, a, w, gamma;
} pil_kapitza;

/* state (theta, omega); use with pil_rk4_step (n = 2) */
void pil_kapitza_deriv(double t, const double *y, double *dydt, void *user);

double pil_kapitza_pivot(const pil_kapitza *p, double t); /* a cos(w t) */
double pil_kapitza_veff(const pil_kapitza *p, double theta);
double pil_kapitza_critical_amp(double g, double l,
                                double w); /* sqrt(2 g l) / w */
/* Small-oscillation frequency of slow motion about theta = 0
 * sqrt(a^2 w^2 / (2 l^2) - g / l), or 0 when inverted position is unstable */
double pil_kapitza_slow_freq(const pil_kapitza *p);

#ifdef __cplusplus
}
#endif

#endif
