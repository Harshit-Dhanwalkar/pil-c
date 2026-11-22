/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_spherical.h - spherical pendulum as a Lagrangian system (use with
 * pil_add_lagrangian, two generalised coordinates q = (theta, phi)) */
#ifndef pil_SPHERICAL_H
#define pil_SPHERICAL_H

#ifdef __cplusplus
extern "C" {
#endif

/* Bob of mass m on a rigid rod of length l; theta is measured from downward
 * vertical, phi is azimuth:
 *   T = 1/2 m l^2 (theta'^2 + sin^2(theta) phi'^2),   V = -m g l cos(theta) */
typedef struct {
  double m, l, g;
} pil_spherical;

double pil_spherical_T(const double *q, const double *qd, void *user);
double pil_spherical_V(const double *q, void *user);

/* Bob position with z up and pivot at origin */
void pil_spherical_bob(const pil_spherical *p, const double *q, double *x,
                       double *y, double *z);
/* Conserved angular momentum about vertical, m l^2 sin^2(theta) phi'
 * (phi is cyclic) */
double pil_spherical_lz(const pil_spherical *p, const double *q,
                        const double *qd);
/* phi' that makes circular (conical) pendulum at theta: sqrt(g / (l
 * cos(theta))) */
double pil_spherical_conical_rate(const pil_spherical *p, double theta);

#ifdef __cplusplus
}
#endif

#endif
