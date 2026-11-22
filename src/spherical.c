/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

#include "pil/pil_spherical.h"
#include <math.h>

double pil_spherical_T(const double *q, const double *qd, void *user) {
  const pil_spherical *p = user;
  double s = sin(q[0]);

  return 0.5 * p->m * p->l * p->l * (qd[0] * qd[0] + s * s * qd[1] * qd[1]);
}

double pil_spherical_V(const double *q, void *user) {
  const pil_spherical *p = user;

  return -p->m * p->g * p->l * cos(q[0]);
}

void pil_spherical_bob(const pil_spherical *p, const double *q, double *x,
                       double *y, double *z) {
  *x = p->l * sin(q[0]) * cos(q[1]);
  *y = p->l * sin(q[0]) * sin(q[1]);
  *z = -p->l * cos(q[0]);
}

double pil_spherical_lz(const pil_spherical *p, const double *q,
                        const double *qd) {
  double s = sin(q[0]);

  return p->m * p->l * p->l * s * s * qd[1];
}

double pil_spherical_conical_rate(const pil_spherical *p, double theta) {
  return sqrt(p->g / (p->l * cos(theta)));
}
