/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_particle.h - charged particles in prescribed electromagnetic fields */
#ifndef pil_PARTICLE_H
#define pil_PARTICLE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Charged particle in 3D with given fields (any E(x, t), B(x, t)), advanced
 * with Boris pusher: second order and |v| is conserved exactly in a pure
 * magnetic field
 * q_m is charge over mass */
typedef void (*pil_field3_fn)(const double x[3], double t, double B[3],
                              double E[3], void *user);
typedef struct {
  double x[3], v[3], q_m;
} pil_particle3;
void pil_particle3_step(pil_particle3 *p, pil_field3_fn field, void *user,
                        double t, double dt);

/* First adiabatic invariant mu = v_perp^2 / (2 |B|) of particle in field at its
 * position (conserved while B changes slowly over a gyration) */
double pil_particle3_moment(const pil_particle3 *p, pil_field3_fn field,
                            void *user, double t);

/* Ready-made fields (pass parameter struct as `user`): uniform B along z;
 * magnetic bottle B_z = b0 (1 + z^2/L^2), B_rho = -b0 rho z / L^2 (two
 * Helmholtz-like coils, field strongest at ends, div B = 0 to first order);
 * dipole of equatorial strength b0 at radius r0 about z axis, B = b0 r0^3 / r^5
 * (-3 z x, -3 z y, r^2 - 3 z^2) */
typedef struct {
  double b0;
} pil_uniform_b;
typedef struct {
  double b0;
  double length;
} pil_bottle;
typedef struct {
  double b0;
  double r0;
} pil_dipole;
void pil_field_uniform_b(const double x[3], double t, double B[3], double E[3],
                         void *user);
void pil_field_bottle(const double x[3], double t, double B[3], double E[3],
                      void *user);
void pil_field_dipole(const double x[3], double t, double B[3], double E[3],
                      void *user);

#ifdef __cplusplus
}
#endif

#endif
