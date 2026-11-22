/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

#include "pil/pil_rotframe.h"
#include <math.h>

void pil_foucault_deriv(double t, const double *y, double *dydt, void *user) {
  (void)t;
  const pil_foucault *f = user;
  dydt[0] = y[2];
  dydt[1] = y[3];
  dydt[2] = -f->w0 * f->w0 * y[0] + 2.0 * f->omega_z * y[3];
  dydt[3] = -f->w0 * f->w0 * y[1] - 2.0 * f->omega_z * y[2];
}

double pil_plane_angle(const double *x, const double *y, int i0, int i1) {
  double sxx = 0.0;
  double syy = 0.0;
  double sxy = 0.0;
  for (int j = i0; j <= i1; j++) {
    sxx += x[j] * x[j];
    syy += y[j] * y[j];
    sxy += x[j] * y[j];
  }

  return 0.5 * atan2(2.0 * sxy, sxx - syy);
}
