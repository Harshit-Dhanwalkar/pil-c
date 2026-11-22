/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

#include "pil/pil_orbit.h"
#include <math.h>

void pil_kepler_step(double s[4], double mu, double dt) {
  double r = hypot(s[0], s[1]);
  double r3 = r * r * r;
  s[2] -= 0.5 * dt * mu * s[0] / r3;
  s[3] -= 0.5 * dt * mu * s[1] / r3;
  s[0] += dt * s[2];
  s[1] += dt * s[3];
  r = hypot(s[0], s[1]);
  r3 = r * r * r;
  s[2] -= 0.5 * dt * mu * s[0] / r3;
  s[3] -= 0.5 * dt * mu * s[1] / r3;
}

double pil_kepler_energy(const double s[4], double mu) {
  return 0.5 * (s[2] * s[2] + s[3] * s[3]) - mu / hypot(s[0], s[1]);
}

double pil_kepler_angmom(const double s[4]) {
  return s[0] * s[3] - s[1] * s[2];
}

void pil_hohmann(double mu, double r1, double r2, double *dv1, double *dv2,
                 double *tof) {
  double a = 0.5 * (r1 + r2);
  if (dv1)
    *dv1 = sqrt(mu / r1) * (sqrt(2.0 * r2 / (r1 + r2)) - 1.0);
  if (dv2)
    *dv2 = sqrt(mu / r2) * (1.0 - sqrt(2.0 * r1 / (r1 + r2)));
  if (tof)
    *tof = M_PI * sqrt(a * a * a / mu);
}

void pil_burn_prograde(double s[4], double dv) {
  double sp = hypot(s[2], s[3]);
  s[2] += dv * s[2] / sp;
  s[3] += dv * s[3] / sp;
}
