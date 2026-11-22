/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_orbit.h - two-body orbits: Kepler motion, Hohmann transfers */
#ifndef pil_ORBIT_H
#define pil_ORBIT_H

#ifdef __cplusplus
extern "C" {
#endif

/* Two-body motion about a point mass with gravitational parameter mu, state
 * (x, y, vx, vy) in orbital plane. Kick-drift-kick leapfrog: symplectic, so
 * energy error stays bounded */
void pil_kepler_step(double s[4], double mu, double dt);

/* Specific orbital energy v^2/2 - mu/r and specific angular momentum
 * x vy - y vx */
double pil_kepler_energy(const double s[4], double mu);
double pil_kepler_angmom(const double s[4]);

/* Hohmann transfer from a circular orbit r1 to a circular orbit r2 (either
 * way round): two tangential velocity changes (negative when braking) and
 * transfer time, half period of ellipse with a = (r1 + r2) / 2 */
void pil_hohmann(double mu, double r1, double r2, double *dv1, double *dv2,
                 double *tof);

/* Apply a velocity change dv along current velocity (tangential burn) */
void pil_burn_prograde(double s[4], double dv);

#ifdef __cplusplus
}
#endif

#endif
