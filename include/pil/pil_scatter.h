/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_scatter.h - Rutherford (Coulomb) scattering of a particle off a fixed
 * repulsive centre. Force per unit mass K / r^2 away from origin (K = k q1 q2 /
 * m); a particle arriving with speed v0 and impact parameter b is deflected
 * through
 *   tan(theta / 2) = kappa / b,     kappa = K / v0^2,
 * and comes no closer than r_min = kappa + sqrt(kappa^2 + b^2) */
#ifndef pil_SCATTER_H
#define pil_SCATTER_H

#ifdef __cplusplus
extern "C" {
#endif

/* Analytic deflection angle in [0, pi] and distance of closest approach */
double pil_scatter_angle(double K, double v0, double b);
double pil_scatter_rmin(double K, double v0, double b);
/* Rutherford differential cross section per solid angle for a beam in 3D,
 * (b / sin(theta)) |db / dtheta| = (kappa / 2)^2 / sin^4(theta / 2) and
 * impact parameter that gives deflection theta */
double pil_scatter_dsigma(double K, double v0, double theta);
double pil_scatter_impact(double K, double v0, double theta);

/* Follow a particle that comes in from x0 < 0 along +x as if from infinity with
 * speed v0 and impact parameter b (start speed and height are corrected for
 * potential already felt at distance |x0|), with leapfrog steps of dt until it
 * is |x0| away from centre again (or nmax steps)
 * Fills path[i] = (x, y) (capacity nmax + 1), returns number of points;
 * measured deflection angle (direction of final velocity, positive when pushed
 * towards +y) goes to *theta and least distance to *rmin (either may be NULL)
 */
int pil_scatter_trace(double K, double v0, double b, double x0, double dt,
                      int nmax, double (*path)[2], double *theta, double *rmin);

#ifdef __cplusplus
}
#endif

#endif
