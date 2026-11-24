/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

#include "pil/pil_optics.h"
#include <math.h>

int pil_snell(double n1, double n2, double ti, double *tt) {
  double s = n1 * sin(ti) / n2;
  if (fabs(s) > 1.0) {
    return 0;
  }
  *tt = asin(s);

  return 1;
}

double pil_critical_angle(double n1, double n2) {
  return n1 > n2 ? asin(n2 / n1) : NAN;
}

double pil_brewster_angle(double n1, double n2) { return atan(n2 / n1); }

void pil_reflect(const double d[2], const double nrm[2], double out[2]) {
  double k = 2.0 * (d[0] * nrm[0] + d[1] * nrm[1]);
  out[0] = d[0] - k * nrm[0];
  out[1] = d[1] - k * nrm[1];
}

int pil_refract(const double d[2], const double nrm[2], double n1, double n2,
                double out[2]) {
  double eta = n1 / n2;
  double c = -(d[0] * nrm[0] + d[1] * nrm[1]);
  double k = 1.0 - eta * eta * (1.0 - c * c);
  if (k < 0.0) {
    pil_reflect(d, nrm, out);

    return 0;
  }

  double f = eta * c - sqrt(k);
  out[0] = eta * d[0] + f * nrm[0];
  out[1] = eta * d[1] + f * nrm[1];

  return 1;
}

double pil_fresnel_rs(double n1, double n2, double ti) {
  double tt;
  if (!pil_snell(n1, n2, ti, &tt)) {
    return 1.0;
  }
  double a = n1 * cos(ti);
  double b = n2 * cos(tt);
  double r = (a - b) / (a + b);

  return r * r;
}

double pil_fresnel_rp(double n1, double n2, double ti) {
  double tt;
  if (!pil_snell(n1, n2, ti, &tt)) {
    return 1.0;
  }
  double a = n2 * cos(ti);
  double b = n1 * cos(tt);
  double r = (a - b) / (a + b);

  return r * r;
}

double pil_prism_deviation(double apex, double n, double ti) {
  double t1;
  if (!pil_snell(1.0, n, ti, &t1)) {
    return NAN;
  }
  double t2 = apex - t1;
  double t3;
  if (!pil_snell(n, 1.0, t2, &t3)) {
    return NAN;
  }

  return ti + t3 - apex;
}

double pil_prism_min_deviation(double apex, double n) {
  return 2.0 * asin(n * sin(0.5 * apex)) - apex;
}

double pil_lens_image(double f, double dobj) {
  return 1.0 / (1.0 / f - 1.0 / dobj);
}

double pil_lens_magnification(double f, double dobj) {
  return -pil_lens_image(f, dobj) / dobj;
}
