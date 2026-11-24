/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

#include "pil/pil_electrostatic.h"
#include <math.h>

void pil_es_field(const pil_charges *c, double x, double y, double *ex,
                  double *ey) {
  *ex = *ey = 0.0;
  for (int i = 0; i < c->n; i++) {
    double dx = x - c->x[i][0];
    double dy = y - c->x[i][1];
    double r2 = dx * dx + dy * dy;
    double r3 = r2 * sqrt(r2);
    if (r3 < 1e-18) {
      continue;
    }
    *ex += c->q[i] * dx / r3;
    *ey += c->q[i] * dy / r3;
  }
}

double pil_es_potential(const pil_charges *c, double x, double y) {
  double v = 0.0;
  for (int i = 0; i < c->n; i++) {
    double r = hypot(x - c->x[i][0], y - c->x[i][1]);
    if (r > 1e-12)
      v += c->q[i] / r;
  }

  return v;
}

void pil_es_field2(double x, double y, double *ex, double *ey, void *user) {
  pil_es_field(user, x, y, ex, ey);
}

int pil_es_field_line(const pil_charges *c, double x0, double y0, int dir,
                      int skip, double ds, int nmax, double xmax, double ymax,
                      double rcut, double *px, double *py) {
  double x = x0;
  double y = y0;
  int n = 0;
  px[n] = x;
  py[n++] = y;
  for (int i = 0; i < nmax; i++) {
    double ex;
    double ey;
    pil_es_field(c, x, y, &ex, &ey);
    double em = hypot(ex, ey);
    if (em < 1e-12) {
      break;
    }
    double h = dir * ds;
    double mx = x + 0.5 * h * ex / em;
    double my = y + 0.5 * h * ey / em;
    pil_es_field(c, mx, my, &ex, &ey);
    em = hypot(ex, ey);
    if (em < 1e-12) {
      break;
    }
    x += h * ex / em;
    y += h * ey / em;
    px[n] = x;
    py[n++] = y;
    if (fabs(x) > xmax || fabs(y) > ymax) {
      break;
    }
    int stop = 0;
    for (int k = 0; k < c->n; k++) {
      if (k == skip && i < 30) {
        continue; // Just left the emitting charge
      }
      if (hypot(x - c->x[k][0], y - c->x[k][1]) < rcut) {
        stop = 1;
      }
    }

    if (stop) {
      break;
    }
  }

  return n;
}

double pil_es_flux_sphere(const pil_charges *c, double cx, double R,
                          int n_theta) {
  double sum = 0.0;
  double dth = M_PI / n_theta;
  for (int i = 0; i < n_theta; i++) {
    double th = (i + 0.5) * dth;
    double nx = cos(th);
    double ny = sin(th);
    double ex;
    double ey;
    pil_es_field(c, cx + R * nx, R * ny, &ex, &ey);
    sum += (ex * nx + ey * ny) * sin(th) * dth;
  }

  return 2.0 * M_PI * R * R * sum;
}

double pil_es_energy(const pil_charges *c) {
  double u = 0.0;
  for (int i = 0; i < c->n; i++) {
    for (int j = i + 1; j < c->n; j++) {
      u += c->q[i] * c->q[j] /
           hypot(c->x[i][0] - c->x[j][0], c->x[i][1] - c->x[j][1]);
    }
  }

  return u;
}
