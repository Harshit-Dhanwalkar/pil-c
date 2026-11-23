/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

#include "pil/pil_orbit.h"
#include "pil/pil_scatter.h"
#include <math.h>

double pil_scatter_angle(double K, double v0, double b) {
  return 2.0 * atan2(K / (v0 * v0), b);
}

double pil_scatter_rmin(double K, double v0, double b) {
  double kp = K / (v0 * v0);

  return kp + sqrt(kp * kp + b * b);
}

double pil_scatter_dsigma(double K, double v0, double theta) {
  double kp = K / (v0 * v0);
  double s = sin(0.5 * theta);

  return 0.25 * kp * kp / (s * s * s * s);
}

double pil_scatter_impact(double K, double v0, double theta) {
  return K / (v0 * v0) / tan(0.5 * theta);
}

int pil_scatter_trace(double K, double v0, double b, double x0, double dt,
                      int nmax, double (*path)[2], double *theta,
                      double *rmin) {
  // Start at finite distance with energy and angular momentum of ideal
  // infinite-distance beam: v_s^2 = v0^2 - 2 K / r, y = b v0 / v_s
  double y = b, vs = v0;
  for (int it = 0; it < 6; it++) {
    vs = sqrt(v0 * v0 - 2.0 * K / hypot(x0, y));
    y = b * v0 / vs;
  }

  double s[4] = {x0, y, vs, 0.0};
  double rm = hypot(x0, y);
  int n = 0;
  path[n][0] = s[0];
  path[n++][1] = s[1];
  for (int i = 0; i < nmax; i++) {
    pil_kepler_step(s, -K, dt); // Repulsion = attraction with mu = -K
    path[n][0] = s[0];
    path[n++][1] = s[1];
    double r = hypot(s[0], s[1]);
    if (r < rm) {
      rm = r;
    }
    if (r > fabs(x0) && s[0] * s[2] + s[1] * s[3] > 0.0) {
      break;
    }
  }

  if (theta) {
    *theta = atan2(s[3], s[2]);
  }
  if (rmin) {
    *rmin = rm;
  }

  return n;
}
