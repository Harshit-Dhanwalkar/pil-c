/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_nbody.h - Newtonian gravity of point masses in plane */
#ifndef pil_NBODY_H
#define pil_NBODY_H

#ifdef __cplusplus
extern "C" {
#endif

#define PIL_NBODY_MAX 16

/* n bodies (n <= PIL_NBODY_MAX) with gravitational constant g and optional
 * Plummer softening eps (0 = exact point masses) */
typedef struct {
  int n;
  double g, eps;
  double m[PIL_NBODY_MAX];
  double x[PIL_NBODY_MAX][2], v[PIL_NBODY_MAX][2];
} pil_nbody;

/* One kick-drift-kick leapfrog step: symplectic and time-reversible, so energy
 * error stays bounded */
void pil_nbody_step(pil_nbody *s, double dt);

/* Total energy (kinetic + potential), momentum and angular momentum about
 * origin, centre of mass position */
double pil_nbody_energy(const pil_nbody *s);
void pil_nbody_momentum(const pil_nbody *s, double p[2]);
double pil_nbody_angmom(const pil_nbody *s);
void pil_nbody_com(const pil_nbody *s, double c[2]);

/* Figure-eight choreography of Chenciner and Montgomery (2000): three equal
 * masses (m = 1, g = 1) chase one another along a single figure-eight curve
 * with period PIL_FIGURE_EIGHT_PERIOD */
#define PIL_FIGURE_EIGHT_PERIOD 6.32591398
void pil_nbody_figure_eight(pil_nbody *s);

#ifdef __cplusplus
}
#endif

#endif
