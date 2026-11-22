/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_biot.h - Biot-Savart fields of current loops */
#ifndef pil_BIOT_H
#define pil_BIOT_H

#ifdef __cplusplus
extern "C" {
#endif

/* Field of a circular loop of radius r carrying unit current (mu0 / 4 pi = 1),
 * axis along z with centre at height z0, summed over nseg straight segments
 * (200 is plenty)
 * Evaluated at (x, 0, z) in plane through axis
 * Returns Bx and Bz (By = 0 there) */
void pil_biot_loop(double r, double z0, double x, double z, int nseg,
                   double *bx, double *bz);

/* Two coaxial loops of radius 1 at z = -d/2 and +d/2 (a Helmholtz pair when
 * d = 1) */
void pil_biot_pair(double d, double x, double z, int nseg, double *bx,
                   double *bz);

/* Follow a field line through (x0, z0) for up to nmax steps of length ds in
 * (x, z) plane, in direction of B (dir = +1) or against it (-1)
 * Stops when line leaves |x| < xmax, |z| < zmax or comes back to its start
 * Points go to px/pz (capacity nmax + 1)
 * Returns how many */
typedef void (*pil_field2_fn)(double x, double z, double *bx, double *bz,
                              void *user);
int pil_field_line(pil_field2_fn field, void *user, double x0, double z0,
                   int dir, double ds, int nmax, double xmax, double zmax,
                   double *px, double *pz);

/* Exact on-axis field of one loop of radius r:
* 2 pi r^2 / (r^2 + (z - z0)^2)^1.5 */
double pil_biot_loop_axis(double r, double z0, double z);

#ifdef __cplusplus
}
#endif

#endif
