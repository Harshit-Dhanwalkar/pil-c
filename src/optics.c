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

double pil_lens_ray_slope(double f, double y, double s) { return s - y / f; }

void pil_lens_principal_slopes(double f, double dobj, double h, double s[3]) {
  s[0] = 0.0;
  s[1] = -h / dobj;
  s[2] = h / (f - dobj);
}

double pil_cauchy_index(double a, double b, double lambda_um) {
  return a + b / (lambda_um * lambda_um);
}

void pil_wavelength_rgb(double l, double rgb[3]) {
  double r = 0.0, g = 0.0, b = 0.0;
  if (l >= 380 && l < 440) {
    r = (440 - l) / 60.0;
    b = 1.0;
  } else if (l < 490) {
    g = (l - 440) / 50.0;
    b = 1.0;
  } else if (l < 510) {
    g = 1.0;
    b = (510 - l) / 20.0;
  } else if (l < 580) {
    r = (l - 510) / 70.0;
    g = 1.0;
  } else if (l < 645) {
    r = 1.0;
    g = (645 - l) / 65.0;
  } else if (l <= 780) {
    r = 1.0;
  }
  double f = l < 420 ? 0.3 + 0.7 * (l - 380) / 40.0
                     : (l > 700 ? 0.3 + 0.7 * (780 - l) / 80.0 : 1.0);
  if (f < 0.0) {
    f = 0.0;
  }
  rgb[0] = r * f;
  rgb[1] = g * f;
  rgb[2] = b * f;
}

/* intersection of the ray p + t d (t > eps) with edge a-b; returns t or -1 */
static double ray_edge(const double p[2], const double d[2], const double a[2],
                       const double b[2]) {
  double ex = b[0] - a[0];
  double ey = b[1] - a[1];
  double den = d[0] * ey - d[1] * ex;
  if (fabs(den) < 1e-14) {
    return -1.0;
  }
  double qx = a[0] - p[0];
  double qy = a[1] - p[1];
  double t = (qx * ey - qy * ex) / den;
  double u = (qx * d[1] - qy * d[0]) / den;
  if (t < 1e-9 || u < -1e-9 || u > 1.0 + 1e-9) {
    return -1.0;
  }

  return t;
}

int pil_glass_trace(const double (*v)[2], int nv, double n,
                    const double start[2], const double dir[2], double tail,
                    double (*pts)[2], int max_pts, double out_dir[2]) {
  double p[2] = {start[0], start[1]};
  double d[2] = {dir[0], dir[1]};
  double dl = hypot(d[0], d[1]);
  d[0] /= dl;
  d[1] /= dl;
  int np = 0;
  int inside = 0;
  int last_edge = -1;
  pts[np][0] = p[0];
  pts[np++][1] = p[1];
  for (int it = 0; it < 8 && np < max_pts - 1; it++) {
    double best = 1e30;
    int be = -1;
    for (int e = 0; e < nv; e++) {
      if (e == last_edge) {
        continue;
      }
      double t = ray_edge(p, d, v[e], v[(e + 1) % nv]);
      if (t > 0 && t < best) {
        best = t;
        be = e;
      }
    }

    if (be < 0) {
      break;
    }
    p[0] += best * d[0];
    p[1] += best * d[1];
    pts[np][0] = p[0];
    pts[np++][1] = p[1];
    double ex = v[(be + 1) % nv][0] - v[be][0];
    double ey = v[(be + 1) % nv][1] - v[be][1];
    double el = hypot(ex, ey);
    double out_n[2] = {ey / el, -ex / el}; // Outward for CCW
    double nrm[2];
    double o[2];
    if (!inside) { // Arriving from outside: normal points back out
      nrm[0] = out_n[0];
      nrm[1] = out_n[1];
      pil_refract(d, nrm, 1.0, n, o);
      inside = 1;
    } else {
      nrm[0] = -out_n[0];
      nrm[1] = -out_n[1];
      inside = pil_refract(d, nrm, n, 1.0, o) ? 0 : 1;
    }

    d[0] = o[0];
    d[1] = o[1];
    last_edge = be;
    if (!inside) {
      break;}
  }

  if (!inside && np > 1) {
    pts[np][0] = p[0] + tail * d[0];
    pts[np++][1] = p[1] + tail * d[1];
  } else if (np > 1) {
    np = 0; // Trapped (should not happen for convex glass in a few bounces)
  } else {
    np = 0;
  }

  if (out_dir && np) {
    out_dir[0] = d[0];
    out_dir[1] = d[1];
  }

  return np;
}
