/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_electrostatic.h - static point charges: field, potential, field lines,
 * Gauss's law. Units with E = q / r^2 (Gaussian-like, 4 pi eps0 = 1), so flux
 * of E out of a closed surface is 4 pi Q_enclosed
 * Charges sit in a plane (x, y); fields are 3D fields evaluated in that plane
 */
#ifndef pil_ELECTROSTATIC_H
#define pil_ELECTROSTATIC_H

#ifdef __cplusplus
extern "C" {
#endif

#define PIL_ES_MAX 16

typedef struct {
  int n;
  double q[PIL_ES_MAX];
  double x[PIL_ES_MAX][2];
} pil_charges;

/* Field E = sum q (r - r_i) / |r - r_i|^3 and potential V = sum q / |r - r_i|
 */
void pil_es_field(const pil_charges *c, double x, double y, double *ex,
                  double *ey);
double pil_es_potential(const pil_charges *c, double x, double y);
/* Same field with (x, y, *ex, *ey, user) signature of pil_field2_fn in
 * pil_biot.h, user = const pil_charges * */
void pil_es_field2(double x, double y, double *ex, double *ey, void *user);

/* Electric field line from (x0, y0): follows +E when dir = +1 and -E when dir =
 * -1 in steps of length ds, until it leaves |x| < xmax, |y| < ymax, gets within
 * rcut of a charge other than one it started at, or after nmax steps
 * Points go to px/py (capacity nmax + 1); returns number of points
 * Start index of charge that emits line is `skip` (-1 for none) */
int pil_es_field_line(const pil_charges *c, double x0, double y0, int dir,
                      int skip, double ds, int nmax, double xmax, double ymax,
                      double rcut, double *px, double *py);

/* Equipotential curve through (x0, y0): steps of length ds along tangent
 * (perpendicular to E) with a Newton correction back onto V = V(x0, y0)
 * Stops when it leaves |x| < xmax, |y| < ymax, gets within rcut of a charge,
 * closes on its start (start point is then repeated as last point) or after
 * nmax steps
 * An open curve is followed in both directions from start (points ordered along
 * it)
 * px/py need 2 nmax + 2 points; returns number of points */
int pil_es_equipotential(const pil_charges *c, double x0, double y0, double ds,
                         int nmax, double xmax, double ymax, double rcut,
                         double *px, double *py);

/* Flux of E through a sphere of radius R centred at (cx, 0, 0), for charges
 * lying on x axis (axially symmetric): 2 pi R^2 int E . n sin(theta) dtheta
 * with n_theta midpoint samples. Equals 4 pi times charge inside */
double pil_es_flux_sphere(const pil_charges *c, double cx, double R,
                          int n_theta);

/* Electrostatic energy sum_{i<j} q_i q_j / r_ij */
double pil_es_energy(const pil_charges *c);

#ifdef __cplusplus
}
#endif

#endif
