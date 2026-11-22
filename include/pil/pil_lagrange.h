/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/*
 * pil_lagrange.h - systems written in generalised coordinates, plus a small 3D
 * projection helper for drawing them
 *
 * Give kinetic energy T(q, qdot) and potential V(q) of figure; Euler-Lagrange
 * equations
 *
 *   d/dt (dT/dqdot) - dT/dq + dV/dq = Q
 *
 * are formed numerically (central differences) and integrated with RK4 as world
 * steps, so plots, recording, reset and animation window work as for everything
 * else. For drawing use pil_world_add_overlay() and projection below
 *
 *   static double T(const double *q, const double *qd, void *u) {...}
 *   static double V(const double *q, void *u) {...}
 *   pil_lsys *s = pil_add_lagrangian(w, 2, T, V, &params);
 *   pil_lsys_set_state(s, (double[]){0.2, 0.3}, (double[]){0, 0});
 *
 * T may contain terms linear in qdot (a rotating frame): they produce Coriolis
 * / gyroscopic terms
 * Up to PIL_LSYS_MAX_DOF coordinates
 */

#ifndef PIL_LAGRANGE_H
#define PIL_LAGRANGE_H

#include "pil_core.h"
#include "pil_particle.h" /* pil_particle3, pil_field3_fn, pil_particle3_step */

#ifdef __cplusplus
extern "C" {
#endif

#define PIL_LSYS_MAX_DOF 12

typedef struct pil_lsys pil_lsys;

typedef double (*pil_kinetic_fn)(const double *q, const double *qd, void *user);
typedef double (*pil_potential_fn)(const double *q, void *user);
/* Generalised forces Q (length n), may depend on time (drives, damping) */
typedef void (*pil_genforce_fn)(const double *q, const double *qd, double t,
                                double *Q, void *user);

pil_lsys *pil_add_lagrangian(pil_world *w, int n, pil_kinetic_fn T,
                             pil_potential_fn V, void *user);
void pil_lsys_set_forces(pil_lsys *s, pil_genforce_fn Q);
/* Initial state; also state pil_world_reset() returns to */
void pil_lsys_set_state(pil_lsys *s, const double *q, const double *qd);
const double *pil_lsys_q(const pil_lsys *s);
const double *pil_lsys_qd(const pil_lsys *s);
int pil_lsys_dof(const pil_lsys *s);
/* Mass matrix M = d2T/dqdot2 (row-major n x n) at current state */
void pil_lsys_mass_matrix(const pil_lsys *s, double *M);
/* Jacobi integral h = qdot . dT/dqdot - T + V: energy when T is quadratic, and
 * still conserved in a steadily rotating frame */
double pil_lsys_jacobi(const pil_lsys *s);

/* Free rigid body (no torque) with principal moments I[0..2] about body axes
 * 1, 2, 3, in Euler angles q = (phi, theta, psi), z-x-z: R = Rz(phi) Rx(theta)
 * Rz(psi), whose columns are body axes in space. Pass pil_rigid3_T as T  (and
 * no V) with a pil_rigid3 as `user`
 * Z axis of space is where angular momentum is put by pil_rigid3_state, so
 * theta is angle from L to axis 3 and stays clear of its 0 / pi singularity
 * unless body spins exactly about axis 3 (or about -3) */
typedef struct {
  double I[3];
} pil_rigid3;
double pil_rigid3_T(const double *q, const double *qd, void *user);
/* R[row][col]: column i = body axis e_i in space */
void pil_rigid3_axes(const double *q, double R[3][3]);
/* Angular velocity, body-frame components */
void pil_rigid3_omega(const double *q, const double *qd, double w_body[3]);
/* Angles q and rates qd for angular velocity w_body (body frame) when angular
 * momentum I w is along space z; phi0 (rotation about L) is free */
void pil_rigid3_state(const pil_rigid3 *p, const double w_body[3], double phi0,
                      double *q, double *qd);

/* Orthographic view of a z-up 3D scene: projected point is in same world units
 * as rest of figure, so annotations can use it. yaw turns scene about z, pitch
 * tilts it towards viewer */
typedef struct {
  double ox, oy; /* where 3D origin lands */
  double yaw, pitch, scale;
} pil_proj3;

pil_proj3 pil_proj3_make(double ox, double oy, double yaw, double pitch,
                         double scale);
pil_vec2 pil_proj3_point(const pil_proj3 *p, double x, double y, double z);
/* Depth along viewing direction (larger = nearer), to order drawing */
double pil_proj3_depth(const pil_proj3 *p, double x, double y, double z);

#ifdef __cplusplus
}
#endif
#endif
