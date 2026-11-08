/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* lagrange.c - Euler-Lagrange equations from T(q, qdot) and V(q)
 *
 * With p_i = dT/dqdot_i equations read
 *   M qddot = Q - dV/dq + dT/dq - (dp/dq) qdot,  M_ij = dp_i/dqdot_j
 * Every derivative is a central difference of T or V, so user writes only
 * energies
 * Step is RK4 on (q, qdot)
 */

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

// rotate about z by yaw, then tilt z axis towards viewer by pitch:
//  screen x = x',
//  screen y = z \cos(pitch) + y' sin(pitch)
//  depth = y' cos - z sin(pitch)
// with y' pointing away from viewer
pil_vec2 pil_proj3_point(const pil_proj3 *p, double x, double y, double z) {
  double cy = cos(p->yaw), sy = sin(p->yaw);
  double xr = cy * x - sy * y, yr = sy * x + cy * y;

  return pil_v2(p->ox + p->scale * xr,
                p->oy + p->scale * (z * cos(p->pitch) + yr * sin(p->pitch)));
}

double pil_proj3_depth(const pil_proj3 *p, double x, double y, double z) {
  double cy = cos(p->yaw), sy = sin(p->yaw);
  double yr = sy * x + cy * y;

  return yr * cos(p->pitch) - z * sin(p->pitch);
}
