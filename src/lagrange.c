/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* lagrange.c - Euler-Lagrange equations from T(q, qdot) and V(q)
 *
 * With p_i = dT/dqdot_i equations read
 *   M qddot = Q - dV/dq + dT/dq - (dp/dq) qdot
 *   M_ij    = dp_i/dqdot_j
 * Every derivative is a central difference of T or V, so user writes only
 * energies
 * Step is RK4 on (q, qdot) */

#include "pil/pil_lagrange.h"
#include "pil_internal.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAXN PIL_LSYS_MAX_DOF

struct pil_lsys {
  pil_world *world;
  int n;
  pil_kinetic_fn T;
  pil_potential_fn V;
  pil_genforce_fn Q;
  void *user;
  double q[MAXN], qd[MAXN], q0[MAXN], qd0[MAXN];
  double max_h;
};

static double step_for(double x, double rel) { return rel * (1.0 + fabs(x)); }

// p = dT/dqdot
static void momenta(const pil_lsys *s, const double *q, const double *qd,
                    double *p) {
  double v[MAXN];
  memcpy(v, qd, (size_t)s->n * sizeof *v);
  for (int i = 0; i < s->n; i++) {
    double h = step_for(qd[i], 1e-4);
    v[i] = qd[i] + h;
    double tp = s->T(q, v, s->user);
    v[i] = qd[i] - h;
    double tm = s->T(q, v, s->user);
    v[i] = qd[i];
    p[i] = (tp - tm) / (2 * h);
  }
}

void pil_lsys_mass_matrix(const pil_lsys *s, double *M) {
  int n = s->n;
  double v[MAXN], pp[MAXN], pm[MAXN];
  memcpy(v, s->qd, (size_t)n * sizeof *v);
  for (int j = 0; j < n; j++) {
    double h = step_for(s->qd[j], 1e-3);
    v[j] = s->qd[j] + h;
    momenta(s, s->q, v, pp);
    v[j] = s->qd[j] - h;
    momenta(s, s->q, v, pm);
    v[j] = s->qd[j];
    for (int i = 0; i < n; i++) {
      M[i * n + j] = (pp[i] - pm[i]) / (2 * h);
    }
  }
}

// qdd for state (q, qd) at time t
static int accel(const pil_lsys *s, const double *q, const double *qd, double t,
                 double *qdd) {
  int n = s->n;
  pil_lsys tmp = *s;
  memcpy(tmp.q, q, (size_t)n * sizeof *q);
  memcpy(tmp.qd, qd, (size_t)n * sizeof *qd);

  double M[MAXN * MAXN], rhs[MAXN], Q[MAXN] = {0}, x[MAXN];
  pil_lsys_mass_matrix(&tmp, M);
  if (s->Q) {
    s->Q(q, qd, t, Q, s->user);
  }

  // dp/dq . qdot: derivative of momenta along motion
  double pp[MAXN], pm[MAXN], qa[MAXN], qb[MAXN];
  double sc = 1e-6;
  for (int i = 0; i < n; i++) {
    qa[i] = q[i] + sc * qd[i];
    qb[i] = q[i] - sc * qd[i];
  }
  momenta(s, qa, qd, pp);
  momenta(s, qb, qd, pm);

  for (int i = 0; i < n; i++) {
    double h = step_for(q[i], 1e-6), w[MAXN];
    memcpy(w, q, (size_t)n * sizeof *w);
    w[i] = q[i] + h;
    double tp = s->T(w, qd, s->user), vp = s->V ? s->V(w, s->user) : 0;
    w[i] = q[i] - h;
    double tm = s->T(w, qd, s->user), vm = s->V ? s->V(w, s->user) : 0;
    double dT = (tp - tm) / (2 * h), dV = (vp - vm) / (2 * h);
    rhs[i] = Q[i] - dV + dT - (pp[i] - pm[i]) / (2 * sc);
  }

  // solve M x = rhs by Gaussian elimination with partial pivoting
  for (int i = 0; i < n; i++) {
    x[i] = rhs[i];
  }
  for (int c = 0; c < n; c++) {
    int piv = c;
    for (int r = c + 1; r < n; r++) {
      if (fabs(M[r * n + c]) > fabs(M[piv * n + c])) {
        piv = r;
      }
    }
    if (fabs(M[piv * n + c]) < 1e-14) {
      return 0;
    }
    if (piv != c) {
      for (int k = 0; k < n; k++) {
        double t_ = M[c * n + k];
        M[c * n + k] = M[piv * n + k];
        M[piv * n + k] = t_;
      }
      double t_ = x[c];
      x[c] = x[piv];
      x[piv] = t_;
    }
    for (int r = c + 1; r < n; r++) {
      double f = M[r * n + c] / M[c * n + c];
      for (int k = c; k < n; k++) {
        M[r * n + k] -= f * M[c * n + k];
      }
      x[r] -= f * x[c];
    }
  }
  for (int r = n - 1; r >= 0; r--) {
    for (int k = r + 1; k < n; k++) {
      x[r] -= M[r * n + k] * x[k];
    }
    x[r] /= M[r * n + r];
  }

  memcpy(qdd, x, (size_t)n * sizeof *x);

  return 1;
}

static void lsys_step(pil_lsys *s, double t, double dt) {
  int n = s->n;
  double k1q[MAXN], k1v[MAXN], k2q[MAXN], k2v[MAXN], k3q[MAXN], k3v[MAXN],
      k4q[MAXN], k4v[MAXN], q[MAXN], v[MAXN];

  memcpy(k1q, s->qd, (size_t)n * sizeof *k1q);
  if (!accel(s, s->q, s->qd, t, k1v)) {
    return;
  }

  for (int i = 0; i < n; i++) {
    q[i] = s->q[i] + 0.5 * dt * k1q[i];
    v[i] = s->qd[i] + 0.5 * dt * k1v[i];
  }
  memcpy(k2q, v, (size_t)n * sizeof *k2q);
  accel(s, q, v, t + 0.5 * dt, k2v);
  for (int i = 0; i < n; i++) {
    q[i] = s->q[i] + 0.5 * dt * k2q[i];
    v[i] = s->qd[i] + 0.5 * dt * k2v[i];
  }
  memcpy(k3q, v, (size_t)n * sizeof *k3q);
  accel(s, q, v, t + 0.5 * dt, k3v);
  for (int i = 0; i < n; i++) {
    q[i] = s->q[i] + dt * k3q[i];
    v[i] = s->qd[i] + dt * k3v[i];
  }
  memcpy(k4q, v, (size_t)n * sizeof *k4q);
  accel(s, q, v, t + dt, k4v);

  for (int i = 0; i < n; i++) {
    s->q[i] += dt / 6 * (k1q[i] + 2 * k2q[i] + 2 * k3q[i] + k4q[i]);
    s->qd[i] += dt / 6 * (k1v[i] + 2 * k2v[i] + 2 * k3v[i] + k4v[i]);
  }
}

pil_lsys *pil_add_lagrangian(pil_world *w, int n, pil_kinetic_fn T,
                             pil_potential_fn V, void *user) {
  if (n < 1 || n > MAXN || !T) {
    return NULL;
  }

  pil_lsys *s = calloc(1, sizeof *s);
  if (!s) {
    abort();
  }

  s->world = w;
  s->n = n;
  s->T = T;
  s->V = V;
  s->user = user;
  s->max_h = 1e-3;
  if (w->nlsys + 1 > w->lscap) {
    w->lscap = (w->nlsys + 1) * 2;
    w->lsys = realloc(w->lsys, (size_t)w->lscap * sizeof *w->lsys);
    if (!w->lsys) {
      abort();
    }
  }

  w->lsys[w->nlsys++] = s;

  return s;
}

void pil_lsys_set_forces(pil_lsys *s, pil_genforce_fn Q) { s->Q = Q; }

void pil_lsys_set_state(pil_lsys *s, const double *q, const double *qd) {
  for (int i = 0; i < s->n; i++) {
    s->q[i] = s->q0[i] = q ? q[i] : 0.0;
    s->qd[i] = s->qd0[i] = qd ? qd[i] : 0.0;
  }
}

const double *pil_lsys_q(const pil_lsys *s) { return s->q; }
const double *pil_lsys_qd(const pil_lsys *s) { return s->qd; }
int pil_lsys_dof(const pil_lsys *s) { return s->n; }

double pil_lsys_jacobi(const pil_lsys *s) {
  double p[MAXN], h = 0;
  momenta(s, s->q, s->qd, p);
  for (int i = 0; i < s->n; i++) {
    h += s->qd[i] * p[i];
  }

  return h - s->T(s->q, s->qd, s->user) + (s->V ? s->V(s->q, s->user) : 0.0);
}

// World hooks
void pil_lsys_step_all(pil_world *w, double t, double dt) {
  for (int i = 0; i < w->nlsys; i++) {
    pil_lsys *s = w->lsys[i];
    int sub = (int)ceil(dt / s->max_h);
    if (sub < 1) {
      sub = 1;
    }

    for (int k = 0; k < sub; k++) {
      lsys_step(s, t + k * dt / sub, dt / sub);
    }
  }
}

void pil_lsys_reset_all(pil_world *w) {
  for (int i = 0; i < w->nlsys; i++) {
    pil_lsys *s = w->lsys[i];
    memcpy(s->q, s->q0, sizeof s->q);
    memcpy(s->qd, s->qd0, sizeof s->qd);
  }
}

void pil_lsys_free_all(pil_world *w) {
  for (int i = 0; i < w->nlsys; i++) {
    free(w->lsys[i]);
  }

  free(w->lsys);
  w->lsys = NULL;
  w->nlsys = 0;
}

// 3D view

pil_proj3 pil_proj3_make(double ox, double oy, double yaw, double pitch,
                         double scale) {
  pil_proj3 p = {ox, oy, yaw, pitch, scale};

  return p;
}

// Rotate about z by yaw, then tilt z axis towards viewer by pitch:
//  screen x = x',
//  screen y = z \cos(pitch) + y' sin(pitch)
//  depth = y' cos - z sin(pitch)
// with y' pointing away from viewer
pil_vec2 pil_proj3_point(const pil_proj3 *p, double x, double y, double z) {
  double cy = cos(p->yaw);
  double sy = sin(p->yaw);
  double xr = cy * x - sy * y;
  double yr = sy * x + cy * y;

  return pil_v2(p->ox + p->scale * xr,
                p->oy + p->scale * (z * cos(p->pitch) + yr * sin(p->pitch)));
}

double pil_proj3_depth(const pil_proj3 *p, double x, double y, double z) {
  double cy = cos(p->yaw);
  double sy = sin(p->yaw);
  double yr = sy * x + cy * y;

  return yr * cos(p->pitch) - z * sin(p->pitch);
}

// Free rigid body in z-x-z Euler angles
void pil_rigid3_omega(const double *q, const double *qd, double w[3]) {
  double st = sin(q[1]);
  double ct = cos(q[1]);
  double sp = sin(q[2]);
  double cp = cos(q[2]);
  w[0] = qd[0] * st * sp + qd[1] * cp;
  w[1] = qd[0] * st * cp - qd[1] * sp;
  w[2] = qd[0] * ct + qd[2];
}

double pil_rigid3_T(const double *q, const double *qd, void *user) {
  const pil_rigid3 *p = user;
  double w[3];
  pil_rigid3_omega(q, qd, w);

  return 0.5 * (p->I[0] * w[0] * w[0] + p->I[1] * w[1] * w[1] +
                p->I[2] * w[2] * w[2]);
}

void pil_rigid3_axes(const double *q, double R[3][3]) {
  double cf = cos(q[0]);
  double sf = sin(q[0]);
  double ct = cos(q[1]);
  double st = sin(q[1]);
  double cp = cos(q[2]);
  double sp = sin(q[2]);
  R[0][0] = cf * cp - sf * ct * sp;
  R[0][1] = -cf * sp - sf * ct * cp;
  R[0][2] = sf * st;
  R[1][0] = sf * cp + cf * ct * sp;
  R[1][1] = -sf * sp + cf * ct * cp;
  R[1][2] = -cf * st;
  R[2][0] = st * sp;
  R[2][1] = st * cp;
  R[2][2] = ct;
}

void pil_rigid3_state(const pil_rigid3 *p, const double w[3], double phi0,
                      double *q, double *qd) {
  double L[3] = {p->I[0] * w[0], p->I[1] * w[1], p->I[2] * w[2]};
  double Lm = sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
  // Space z seen from the body is (sin th sin psi, sin th cos psi, cos th)
  double th = acos(L[2] / Lm);
  double ps = atan2(L[0], L[1]);
  double st = sin(th);
  double ct = cos(th);
  double sp = sin(ps);
  double cp = cos(ps);
  double phid = (w[0] * sp + w[1] * cp) / st;
  q[0] = phi0;
  q[1] = th;
  q[2] = ps;
  qd[0] = phid;
  qd[1] = w[0] * cp - w[1] * sp;
  qd[2] = w[2] - phid * ct;
}

void pil_particle3_step(pil_particle3 *p, pil_field3_fn field, void *user,
                        double t, double dt) {
  double B[3] = {0, 0, 0};
  double E[3] = {0, 0, 0};
  field(p->x, t, B, E, user);

  double h = 0.5 * p->q_m * dt;
  double vm[3];
  double tv[3];
  double vp[3];
  for (int i = 0; i < 3; i++) {
    vm[i] = p->v[i] + h * E[i];
    tv[i] = h * B[i];
  }

  double t2 = tv[0] * tv[0] + tv[1] * tv[1] + tv[2] * tv[2];
  double f = 2.0 / (1.0 + t2);
  vp[0] = vm[0] + vm[1] * tv[2] - vm[2] * tv[1];
  vp[1] = vm[1] + vm[2] * tv[0] - vm[0] * tv[2];
  vp[2] = vm[2] + vm[0] * tv[1] - vm[1] * tv[0];
  double sv[3] = {f * tv[0], f * tv[1], f * tv[2]};
  double vq[3] = {vm[0] + vp[1] * sv[2] - vp[2] * sv[1],
                  vm[1] + vp[2] * sv[0] - vp[0] * sv[2],
                  vm[2] + vp[0] * sv[1] - vp[1] * sv[0]};
  for (int i = 0; i < 3; i++) {
    p->v[i] = vq[i] + h * E[i];
    p->x[i] += p->v[i] * dt;
  }
}
