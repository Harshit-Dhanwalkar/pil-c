/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

#include "pil/pil_particle.h"
#include "pil/pil_lagrange.h"
#include <math.h>

// void pil_particle3_step(pil_particle3 *p, pil_field3_fn field, void *user,
//                         double t, double dt) {
//   double B[3] = {0, 0, 0};
//   double E[3] = {0, 0, 0};
//   field(p->x, t, B, E, user);
//
//   double h = 0.5 * p->q_m * dt;
//   double vm[3];
//   double tv[3];
//   double vp[3];
//   for (int i = 0; i < 3; i++) {
//     vm[i] = p->v[i] + h * E[i];
//     tv[i] = h * B[i];
//   }
//
//   double t2 = tv[0] * tv[0] + tv[1] * tv[1] + tv[2] * tv[2];
//   double f = 2.0 / (1.0 + t2);
//   vp[0] = vm[0] + vm[1] * tv[2] - vm[2] * tv[1];
//   vp[1] = vm[1] + vm[2] * tv[0] - vm[0] * tv[2];
//   vp[2] = vm[2] + vm[0] * tv[1] - vm[1] * tv[0];
//   double sv[3] = {f * tv[0], f * tv[1], f * tv[2]};
//   double vq[3] = {vm[0] + vp[1] * sv[2] - vp[2] * sv[1],
//                   vm[1] + vp[2] * sv[0] - vp[0] * sv[2],
//                   vm[2] + vp[0] * sv[1] - vp[1] * sv[0]};
//   for (int i = 0; i < 3; i++) {
//     p->v[i] = vq[i] + h * E[i];
//     p->x[i] += p->v[i] * dt;
//   }
// }

void pil_field_uniform_b(const double x[3], double t, double B[3], double E[3],
                         void *user) {
  (void)x;
  (void)t;
  const pil_uniform_b *f = user;
  B[0] = B[1] = 0.0;
  B[2] = f->b0;
  E[0] = E[1] = E[2] = 0.0;
}

void pil_field_bottle(const double x[3], double t, double B[3], double E[3],
                      void *user) {
  (void)t;
  const pil_bottle *f = user;
  double l2 = f->length * f->length;
  B[0] = -f->b0 * x[0] * x[2] / l2;
  B[1] = -f->b0 * x[1] * x[2] / l2;
  B[2] = f->b0 * (1.0 + x[2] * x[2] / l2);
  E[0] = E[1] = E[2] = 0.0;
}

void pil_field_dipole(const double x[3], double t, double B[3], double E[3],
                      void *user) {
  (void)t;
  const pil_dipole *f = user;
  double r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  double r5 = r2 * r2 * sqrt(r2);
  double k = f->b0 * f->r0 * f->r0 * f->r0 / r5;
  B[0] = -3.0 * k * x[2] * x[0];
  B[1] = -3.0 * k * x[2] * x[1];
  B[2] = k * (r2 - 3.0 * x[2] * x[2]);
  E[0] = E[1] = E[2] = 0.0;
}

double pil_particle3_moment(const pil_particle3 *p, pil_field3_fn field,
                            void *user, double t) {
  double B[3];
  double E[3];
  field(p->x, t, B, E, user);
  double bm = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
  double vpar = (p->v[0] * B[0] + p->v[1] * B[1] + p->v[2] * B[2]) / bm;
  double v2 = p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2];

  return 0.5 * (v2 - vpar * vpar) / bm;
}
