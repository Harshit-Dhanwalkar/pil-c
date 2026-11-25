/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

#include "pil/pil_electrostatic.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

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
        continue; // Just left emitting charge
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

/* One direction (sgn = +1 / -1) along tangent; returns 1 in *closed if curve
 * came back to its start */
static int equi_dir(const pil_charges *c, double x0, double y0, double ds,
                    int nmax, double xmax, double ymax, double rcut, int sgn,
                    double *px, double *py, int *closed) {
  double v0 = pil_es_potential(c, x0, y0), x = x0, y = y0;
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
    double mx = x - sgn * 0.5 * ds * ey / em;
    double my = y + sgn * 0.5 * ds * ex / em;
    pil_es_field(c, mx, my, &ex, &ey);
    em = hypot(ex, ey);
    if (em < 1e-12) {
      break;
    }
    x -= sgn * ds * ey / em;
    y += sgn * ds * ex / em;
    for (int k = 0; k < 2; k++) { // Newton: dx = (V - V0) E / |E|^2
      pil_es_field(c, x, y, &ex, &ey);
      double e2 = ex * ex + ey * ey;
      double dv = pil_es_potential(c, x, y) - v0;
      if (e2 < 1e-24) {
        break;
      }
      x += dv * ex / e2;
      y += dv * ey / e2;
    }

    if (fabs(x) > xmax || fabs(y) > ymax) {
      break;
    }
    int stop = 0;
    for (int k = 0; k < c->n; k++)
      if (hypot(x - c->x[k][0], y - c->x[k][1]) < rcut) {
        stop = 1;
      }
    if (stop) {
      break;
    }
    if (i > 10 && hypot(x - x0, y - y0) < ds) {
      px[n] = x0;
      py[n++] = y0;
      *closed = 1;
      break;
    }

    px[n] = x;
    py[n++] = y;
  }

  return n;
}

int pil_es_equipotential(const pil_charges *c, double x0, double y0, double ds,
                         int nmax, double xmax, double ymax, double rcut,
                         double *px, double *py) {
  int closed = 0;
  int nf = equi_dir(c, x0, y0, ds, nmax, xmax, ymax, rcut, 1, px, py, &closed);
  if (closed) {
    return nf;
  }
  double *bx = malloc(2 * (size_t)(nmax + 2) * sizeof *bx);
  double *by = bx + nmax + 2;
  if (!bx) {
    return nf;
  }
  int unused = 0,
      nb = equi_dir(c, x0, y0, ds, nmax, xmax, ymax, rcut, -1, bx, by, &unused);
  /* Backward part reversed (without shared start point), then forward */
  memmove(px + (nb - 1), px, (size_t)nf * sizeof *px);
  memmove(py + (nb - 1), py, (size_t)nf * sizeof *py);
  for (int i = 1; i < nb; i++) {
    px[nb - 1 - i] = bx[i];
    py[nb - 1 - i] = by[i];
  }

  free(bx);

  return nf + nb - 1;
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
