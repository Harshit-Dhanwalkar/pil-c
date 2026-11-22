/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

#include "pil/pil_biot.h"
#include <math.h>

void pil_biot_loop(double r, double z0, double x, double z, int nseg,
                   double *bx, double *bz) {
  *bx = *bz = 0.0;
  double dph = 2.0 * M_PI / nseg;
  for (int i = 0; i < nseg; i++) {
    double ph = (i + 0.5) * dph;
    double px = r * cos(ph);
    double py = r * sin(ph);
    double dlx = -r * sin(ph) * dph;
    double dly = r * cos(ph) * dph;
    double rx = x - px;
    double ry = -py;
    double rz = z - z0;
    double r3 = pow(rx * rx + ry * ry + rz * rz, 1.5);
    *bx += dly * rz / r3; // (dl x r)_x = dl_y r_z - dl_z r_y
    *bz += (dlx * ry - dly * rx) / r3;
  }
}

void pil_biot_pair(double d, double x, double z, int nseg, double *bx,
                   double *bz) {
  double ax;
  double az;
  double cx;
  double cz;
  pil_biot_loop(1.0, -0.5 * d, x, z, nseg, &ax, &az);
  pil_biot_loop(1.0, 0.5 * d, x, z, nseg, &cx, &cz);
  *bx = ax + cx;
  *bz = az + cz;
}

int pil_field_line(pil_field2_fn field, void *user, double x0, double z0,
                   int dir, double ds, int nmax, double xmax, double zmax,
                   double *px, double *pz) {
  double x = x0;
  double z = z0;
  double bx;
  double bz;
  int n = 0;
  px[n] = x;
  pz[n++] = z;
  for (int i = 0; i < nmax; i++) {
    field(x, z, &bx, &bz, user);
    double bm = hypot(bx, bz);
    if (bm < 1e-12) {
      break;
    }
    double h = dir * ds;
    double mx = x + 0.5 * h * bx / bm;
    double mz = z + 0.5 * h * bz / bm;
    field(mx, mz, &bx, &bz, user);
    bm = hypot(bx, bz);
    x += h * bx / bm;
    z += h * bz / bm;
    if (fabs(x) > xmax || fabs(z) > zmax) {
      break;
    }
    px[n] = x;
    pz[n++] = z;
    if (i > 20 && hypot(x - x0, z - z0) < ds) {
      break;
    }
  }

  return n;
}

double pil_biot_loop_axis(double r, double z0, double z) {
  double d = z - z0;

  return 2.0 * M_PI * r * r / pow(r * r + d * d, 1.5);
}
