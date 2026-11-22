/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

#include "pil/pil_nbody.h"
#include <math.h>

static void accel(const pil_nbody *s, double a[][2]) {
  for (int i = 0; i < s->n; i++) {
    a[i][0] = a[i][1] = 0.0;
  }
  for (int i = 0; i < s->n; i++) {
    for (int j = i + 1; j < s->n; j++) {
      double dx = s->x[j][0] - s->x[i][0];
      double dy = s->x[j][1] - s->x[i][1];
      double r2 = dx * dx + dy * dy + s->eps * s->eps;
      double r3 = r2 * sqrt(r2);
      double fi = s->g * s->m[j] / r3;
      double fj = s->g * s->m[i] / r3;
      a[i][0] += fi * dx;
      a[i][1] += fi * dy;
      a[j][0] -= fj * dx;
      a[j][1] -= fj * dy;
    }
  }
}

void pil_nbody_step(pil_nbody *s, double dt) {
  double a[PIL_NBODY_MAX][2];
  accel(s, a);
  for (int i = 0; i < s->n; i++) {
    s->v[i][0] += 0.5 * dt * a[i][0];
    s->v[i][1] += 0.5 * dt * a[i][1];
    s->x[i][0] += dt * s->v[i][0];
    s->x[i][1] += dt * s->v[i][1];
  }

  accel(s, a);
  for (int i = 0; i < s->n; i++) {
    s->v[i][0] += 0.5 * dt * a[i][0];
    s->v[i][1] += 0.5 * dt * a[i][1];
  }
}

double pil_nbody_energy(const pil_nbody *s) {
  double e = 0.0;
  for (int i = 0; i < s->n; i++) {
    e += 0.5 * s->m[i] * (s->v[i][0] * s->v[i][0] + s->v[i][1] * s->v[i][1]);
    for (int j = i + 1; j < s->n; j++) {
      double dx = s->x[j][0] - s->x[i][0];
      double dy = s->x[j][1] - s->x[i][1];
      e -= s->g * s->m[i] * s->m[j] / sqrt(dx * dx + dy * dy + s->eps * s->eps);
    }
  }

  return e;
}

void pil_nbody_momentum(const pil_nbody *s, double p[2]) {
  p[0] = p[1] = 0.0;
  for (int i = 0; i < s->n; i++) {
    p[0] += s->m[i] * s->v[i][0];
    p[1] += s->m[i] * s->v[i][1];
  }
}

double pil_nbody_angmom(const pil_nbody *s) {
  double l = 0.0;
  for (int i = 0; i < s->n; i++) {
    l += s->m[i] * (s->x[i][0] * s->v[i][1] - s->x[i][1] * s->v[i][0]);
  }

  return l;
}

void pil_nbody_com(const pil_nbody *s, double c[2]) {
  double mt = 0.0;
  c[0] = c[1] = 0.0;
  for (int i = 0; i < s->n; i++) {
    mt += s->m[i];
    c[0] += s->m[i] * s->x[i][0];
    c[1] += s->m[i] * s->x[i][1];
  }

  c[0] /= mt;
  c[1] /= mt;
}

void pil_nbody_figure_eight(pil_nbody *s) {
  s->n = 3;
  s->g = 1.0;
  s->eps = 0.0;
  for (int i = 0; i < 3; i++) {
    s->m[i] = 1.0;
  }
  s->x[0][0] = -0.97000436;
  s->x[0][1] = 0.24308753;
  s->x[1][0] = 0.97000436;
  s->x[1][1] = -0.24308753;
  s->x[2][0] = 0.0;
  s->x[2][1] = 0.0;
  s->v[2][0] = -0.93240737;
  s->v[2][1] = -0.86473146;
  s->v[0][0] = -0.5 * s->v[2][0];
  s->v[0][1] = -0.5 * s->v[2][1];
  s->v[1][0] = s->v[0][0];
  s->v[1][1] = s->v[0][1];
}
