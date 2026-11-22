/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

#include "pil/pil_ode.h"
#include <math.h>

void pil_rk4_step(pil_ode_fn f, int n, double t, double *y, double h,
                  void *user) {
  double k1[16];
  double k2[16];
  double k3[16];
  double k4[16];
  double u[16];
  f(t, y, k1, user);
  for (int i = 0; i < n; i++) {
    u[i] = y[i] + 0.5 * h * k1[i];
  }
  f(t + 0.5 * h, u, k2, user);
  for (int i = 0; i < n; i++) {
    u[i] = y[i] + 0.5 * h * k2[i];
  }
  f(t + 0.5 * h, u, k3, user);
  for (int i = 0; i < n; i++) {
    u[i] = y[i] + h * k3[i];
  }
  f(t + h, u, k4, user);
  for (int i = 0; i < n; i++) {
    y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
  }
}

double pil_wrap_pi(double a) {
  a = fmod(a + M_PI, 2.0 * M_PI);
  if (a < 0.0) {
    a += 2.0 * M_PI;
  }

  return a - M_PI;
}
