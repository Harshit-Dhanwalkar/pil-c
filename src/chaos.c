/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

#include "pil/pil_chaos.h"
#include <math.h>

#include "pil/pil_ode.h"

void pil_driven_pendulum_deriv(double t, const double *y, double *dydt,
                               void *user) {
  const pil_driven_pendulum *p = user;
  dydt[0] = y[1];
  dydt[1] = -y[1] / p->q - sin(y[0]) + p->g * cos(p->w_d * t);
}

void pil_poincare_section(const pil_driven_pendulum *p, double theta0,
                          double omega0, int nper, int steps_per_period,
                          double (*section)[2], double (*path)[2]) {
  double y[2] = {theta0, omega0};
  double period = 2.0 * M_PI / p->w_d;
  double h = period / steps_per_period;
  for (int k = 0; k < nper; k++) {
    for (int i = 0; i < steps_per_period; i++) {
      double t = (k * (double)steps_per_period + i) * h;
      if (i == 0 && section) {
        section[k][0] = pil_wrap_pi(y[0]);
        section[k][1] = y[1];
      }

      if (path) {
        path[k * steps_per_period + i][0] = pil_wrap_pi(y[0]);
        path[k * steps_per_period + i][1] = y[1];
      }

      pil_rk4_step(pil_driven_pendulum_deriv, 2, t, y, h, (void *)p);
    }
  }
}
