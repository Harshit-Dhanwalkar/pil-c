/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

#include "pil/pil_emwave.h"
#include <math.h>

void pil_emwave_field(const pil_emwave *w, double x, double t, double E[3],
                      double B[3]) {
  double ph = w->k * x - w->c * w->k * t + w->phase;
  E[0] = 0.0;
  E[1] = w->e0 * cos(ph);
  E[2] = w->e0 * w->ellipticity * sin(ph);
  B[0] = 0.0;
  B[1] = -E[2] / w->c;
  B[2] = E[1] / w->c;
}

void pil_emwave_poynting(const pil_emwave *w, double x, double t, double S[3]) {
  double E[3];
  double B[3];
  double mu0 = 1.0 / (w->c * w->c);
  pil_emwave_field(w, x, t, E, B);
  S[0] = (E[1] * B[2] - E[2] * B[1]) / mu0;
  S[1] = (E[2] * B[0] - E[0] * B[2]) / mu0;
  S[2] = (E[0] * B[1] - E[1] * B[0]) / mu0;
}

double pil_emwave_energy_density(const pil_emwave *w, double x, double t) {
  double E[3];
  double B[3];
  double mu0 = 1.0 / (w->c * w->c);
  pil_emwave_field(w, x, t, E, B);

  return 0.5 * (E[0] * E[0] + E[1] * E[1] + E[2] * E[2]) +
         0.5 * (B[0] * B[0] + B[1] * B[1] + B[2] * B[2]) / mu0;
}

double pil_emwave_intensity(const pil_emwave *w) {
  return 0.5 * w->c * w->e0 * w->e0 * (1.0 + w->ellipticity * w->ellipticity);
}

double pil_emwave_pressure(const pil_emwave *w, int reflector) {
  return (reflector ? 2.0 : 1.0) * pil_emwave_intensity(w) / w->c;
}
