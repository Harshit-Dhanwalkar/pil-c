/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* solver.c - forces + constraints -> accelerations
 *
 * Equations of motion for maximal coordinates q = (x, y, \theta) per body:
 *      M a + J^T \lambda = F        J a = \gamma
 * Eliminating a gives (small, dense) Schur system
 *      (J M^-1 J^T) \lambda = J M^-1 F - \gamma
 * solved by Gaussian elimination with partial pivoting.
 * Ropes are handled with a tiny active-set loop (a rope may only pull: tension
 * = +\lambda >= 0, force on body = -J^T \lambda)
 */
#include "pil_internal.h"
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define GROW(arr, n, cap)                                                      \
  do {                                                                         \
    if ((n) > (cap)) {                                                         \
      (cap) = (n) * 2;                                                         \
      (arr) = realloc((arr), (size_t)(cap) * sizeof(*(arr)));                  \
      if (!(arr)) {                                                            \
        abort();                                                               \
      }                                                                        \
    }                                                                          \
  } while (0)

static int solve_dense(double *A, double *b, int n) {
  for (int c = 0; c < n; c++) {
    int p = c;
    for (int r = c + 1; r < n; r++) {
      if (fabs(A[r * n + c]) > fabs(A[p * n + c])) {
        p = r;
      }
    }
    if (fabs(A[p * n + c]) < 1e-300) {
      return 0;
    }

    if (p != c) {
      for (int k = 0; k < n; k++) {
        double t = A[c * n + k];
        A[c * n + k] = A[p * n + k];
        A[p * n + k] = t;
      }

      double t = b[c];
      b[c] = b[p];
      b[p] = t;
    }

    for (int r = c + 1; r < n; r++) {
      double f = A[r * n + c] / A[c * n + c];
      if (f == 0.0) {
        continue;
      }
      for (int k = c; k < n; k++) {
        A[r * n + k] -= f * A[c * n + k];
      }
      b[r] -= f * b[c];
    }
  }

  for (int r = n - 1; r >= 0; r--) {
    double s = b[r];
    for (int k = r + 1; k < n; k++) {
      s -= A[r * n + k] * b[k];
    }
    b[r] = s / A[r * n + r];
  }

  return 1;
}

// coupling of two Jacobian blocks that act on same body
static inline double couple(const pil_body *b, const double *j1,
                            const double *j2) {
  return (j1[0] * j2[0] + j1[1] * j2[1]) * b->inv_mass +
         j1[2] * j2[2] * b->inv_inertia;
}

static int build_rows(pil_world *w) {
  int m = 0;
  for (int i = 0; i < w->nels; i++) {
    pil_element *e = w->els[i];
    if (!pil_is_constraint(e)) {
      continue;
    }

    // a constraint between static things carries no force and would only add a
    // singular row
    if (!pil_isdyn(e->a.body) && !pil_isdyn(e->b.body)) {
      continue;
    }
    pil_row r[2];
    int n = pil_constraint_eval(e, w->tcur, r, 0);
    GROW(w->rows, m + n, w->rowcap);
    for (int k = 0; k < n; k++) {
      w->rows[m].el = e;
      w->rows[m].idx = k;
      w->rows[m].r = r[k];
      w->rows[m].lambda = 0;
      m++;
    }
  }

  return m;
}

// K = J M^-1 J^T for rows in w->rows (stored at w->K; m right-hand-side slots
// follow it)
static void build_K(pil_world *w, int m) {
  if (m <= 0) {
    return;
  }

  int need = m * m + m;
  GROW(w->K, need, w->kcap);
  double *K = w->K;
  for (int i = 0; i < m; i++) {
    const pil_rowrec *ri = &w->rows[i];
    pil_body *ai = ri->el->a.body;
    pil_body *bi = ri->el->b.body;

    for (int j = 0; j < m; j++) {
      const pil_rowrec *rj = &w->rows[j];
      pil_body *aj = rj->el->a.body;
      pil_body *bj = rj->el->b.body;
      double v = 0;
      if (pil_isdyn(ai) && ai == aj) {
        v += couple(ai, ri->r.ja, rj->r.ja);
      }
      if (pil_isdyn(ai) && ai == bj) {
        v += couple(ai, ri->r.ja, rj->r.jb);
      }
      if (pil_isdyn(bi) && bi == aj) {
        v += couple(bi, ri->r.jb, rj->r.ja);
      }
      if (pil_isdyn(bi) && bi == bj) {
        v += couple(bi, ri->r.jb, rj->r.jb);
      }

      K[i * m + j] = v;
    }

    K[i * m + i] +=
        1e-10 * (1.0 + K[i * m + i]); // tolerate redundant constraints
  }
}

static void solve_lambda(pil_world *w, int m) {
  if (m <= 0) {
    return;
  }

  build_K(w, m);
  double *K = w->K;
  double *rhs = w->K + (size_t)m * m;
  for (int i = 0; i < m; i++) {
    const pil_rowrec *ri = &w->rows[i];
    pil_body *ai = ri->el->a.body, *bi = ri->el->b.body;
    // rhs = J M^-1 F - \gamma
    double s = -ri->r.gamma;
    if (pil_isdyn(ai)) {
      s += (ri->r.ja[0] * ai->force.x + ri->r.ja[1] * ai->force.y) *
               ai->inv_mass +
           ri->r.ja[2] * ai->torque * ai->inv_inertia;
    }
    if (pil_isdyn(bi)) {
      s += (ri->r.jb[0] * bi->force.x + ri->r.jb[1] * bi->force.y) *
               bi->inv_mass +
           ri->r.jb[2] * bi->torque * bi->inv_inertia;
    }

    rhs[i] = s;
  }

  if (!solve_dense(K, rhs, m)) {
    memset(rhs, 0, (size_t)m * sizeof(double));
  }
  for (int i = 0; i < m; i++) {
    w->rows[i].lambda = rhs[i];
  }
}

// Impacts. A rope that goes taut while its ends are still moving apart must
// stop that motion at once (inelastic), and a body that arrives at a stop must
// rebound with its restitution; otherwise stabilised constraints would let rope
// stretch or body sink before reacting
// Solve for impulses that give impacting rows their post-impact rate (rope:
// zero extension rate, stop: -e times arrival rate) and leave every other
// constraint's rate untouched, then apply them to body velocities
static int stop_hit(const pil_element *e, pil_row *r) {
  if (e->kind != EL_STOP || !pil_isdyn(e->a.body)) {
    return 0;
  }

  pil_stop_row(e, r);

  return r->phi > -1e-9 && r->vel > 1e-6;
}

void pil_project_rope_snaps(pil_world *w) {
  const double vtol = 1e-4;
  int any = 0;
  for (int i = 0; i < w->nels && !any; i++) {
    pil_element *e = w->els[i];
    pil_row r[2];
    if (stop_hit(e, &r[0])) {
      any = 1;
      continue;
    }

    if (!pil_is_rope(e)) {
      continue;
    }

    if (!pil_isdyn(e->a.body) && !pil_isdyn(e->b.body)) {
      continue;
    }

    if (pil_constraint_eval(e, w->tcur, r, 1) > 0 &&
        r[0].phi > -1e-3 * e->length && r[0].vel > vtol) {
      any = 1;
    }
  }

  if (!any) {
    return;
  }

  double T_before = pil_world_kinetic(w);
  int m = 0;
  for (int i = 0; i < w->nels; i++) {
    pil_element *e = w->els[i];
    pil_row r[2];
    if (e->kind == EL_STOP) {
      if (!stop_hit(e, &r[0])) {
        continue;
      }

      GROW(w->rows, m + 1, w->rowcap);
      w->rows[m].el = e;
      w->rows[m].idx = 0;
      w->rows[m].r = r[0];
      w->rows[m].lambda = 0;
      w->rows[m].target = (1.0 + e->k) * r[0].vel;
      m++;
      continue;
    }

    if (!pil_is_constraint(e)) {
      continue;
    }

    if (!pil_isdyn(e->a.body) && !pil_isdyn(e->b.body)) {
      continue;
    }

    int n = pil_constraint_eval(e, w->tcur, r, 1);
    for (int k = 0; k < n; k++) {
      double target = 0;
      if (pil_is_rope(e)) {
        if (!(r[k].phi > -1e-3 * e->length && r[k].vel > vtol)) {
          continue; // slack or not extending: leave alone
        }

        target = r[k].vel;
      }

      GROW(w->rows, m + 1, w->rowcap);
      w->rows[m].el = e;
      w->rows[m].idx = k;
      w->rows[m].r = r[k];
      w->rows[m].lambda = 0;
      w->rows[m].target = target;
      m++;
    }
  }

  for (int it = 0; it < 4; it++) {
    build_K(w, m);
    double *K = w->K;
    double *rhs = w->K + (size_t)m * m;
    for (int i = 0; i < m; i++) {
      rhs[i] = w->rows[i].target;
    }

    if (!solve_dense(K, rhs, m)) {
      return;
    }

    int changed = 0;
    for (int i = 0; i < m; i++) {
      w->rows[i].lambda = rhs[i];
      if (pil_is_rope(w->rows[i].el) && w->rows[i].target != 0 && rhs[i] < 0) {
        w->rows[i].target = 0;
        changed = 1;
      } // a rope cannot push
    }

    if (!changed) {
      break;
    }
  }

  for (int i = 0; i < m; i++) {
    const pil_rowrec *r = &w->rows[i];
    pil_body *a = r->el->a.body, *b = r->el->b.body;
    if (pil_isdyn(a)) {
      a->vel.x -= r->r.ja[0] * r->lambda * a->inv_mass;
      a->vel.y -= r->r.ja[1] * r->lambda * a->inv_mass;
      a->omega -= r->r.ja[2] * r->lambda * a->inv_inertia;
    }

    if (pil_isdyn(b)) {
      b->vel.x -= r->r.jb[0] * r->lambda * b->inv_mass;
      b->vel.y -= r->r.jb[1] * r->lambda * b->inv_mass;
      b->omega -= r->r.jb[2] * r->lambda * b->inv_inertia;
    }
  }

  double lost = T_before - pil_world_kinetic(w);
  if (lost > 0) {
    w->dissipated += lost;
  }
}

void pil_compute_accel(pil_world *w) {
  // 1. applied forces
  for (int i = 0; i < w->nbodies; i++) {
    pil_body *b = w->bodies[i];
    b->force = pil_v2(0, 0);
    b->torque = 0;
    if (b->dof < 0) {
      continue;
    }
    b->force =
        v_add(v_mul(w->g, b->mass), v_mul(b->vel, -w->drag_lin * b->mass));
    b->torque = -w->drag_ang * b->inertia * b->omega;
  }

  for (int i = 0; i < w->nels; i++) {
    if (!pil_is_constraint(w->els[i])) {
      pil_element_apply(w->els[i], w, w->tcur);
    }
  }

  pil_circuits_apply_forces(w); // motors and voice coils

  if (w->grab.b) {
    pil_vec2 p = pil_body_point(w->grab.b, w->grab.local);
    pil_vec2 f = v_sub(v_mul(v_sub(w->grab.target, p), w->grab.k),
                       v_mul(pil_body_velocity_at(w->grab.b, p), w->grab.c));
    pil_apply_force(w->grab.b, p, f);
    // Point damper cannot see rotation about (stationary) grab point, so a
    // dragged body would spin forever; damp it at same rate c/m
    {
      pil_body *gb = w->grab.b;
      pil_vec2 r = v_sub(p, gb->pos);
      pil_apply_torque(gb, -(w->grab.c / gb->mass) *
                               (gb->inertia + gb->mass * v_dot(r, r)) *
                               gb->omega);
    }
  }

  // 2. rope state
  for (int i = 0; i < w->nels; i++) {
    pil_element *e = w->els[i];
    if (!pil_is_rope(e)) {
      continue;
    }

    // RK4 stages can undershoot by ~1e-6; real slackness is caught by \lambda
    // test below
    e->rope_off = pil_rope_length(e) < e->length * (1.0 - 1e-3);
  }

  // 3. constraint multipliers (with rope active-set loop)
  int m = 0;
  for (int it = 0; it < 6; it++) {
    m = build_rows(w);
    if (m == 0) {
      break;
    }

    solve_lambda(w, m);
    int changed = 0;
    for (int i = 0; i < m; i++) {
      pil_element *e = w->rows[i].el;
      if (pil_is_rope(e) && w->rows[i].lambda < -1e-12) {
        e->rope_off = 1;
        changed = 1;
      }
    }

    if (!changed) {
      break;
    }
  }

  // 4. accelerations
  for (int i = 0; i < w->nbodies; i++) {
    pil_body *b = w->bodies[i];
    if (b->dof < 0) {
      b->acc = pil_v2(0, 0);
      b->alpha = 0;
      continue;
    }

    b->acc = v_mul(b->force, b->inv_mass);
    b->alpha = b->torque * b->inv_inertia;
  }

  for (int i = 0; i < w->nels; i++) {
    if (pil_is_constraint(w->els[i])) {
      w->els[i]->last_force = 0;
    }
  }

  for (int i = 0; i < m; i++) {
    pil_rowrec *r = &w->rows[i];
    double lam = r->lambda;
    pil_body *a = r->el->a.body, *b = r->el->b.body;
    if (pil_isdyn(a)) {
      a->acc.x -= r->r.ja[0] * lam * a->inv_mass;
      a->acc.y -= r->r.ja[1] * lam * a->inv_mass;
      a->alpha -= r->r.ja[2] * lam * a->inv_inertia;
    }

    if (pil_isdyn(b)) {
      b->acc.x -= r->r.jb[0] * lam * b->inv_mass;
      b->acc.y -= r->r.jb[1] * lam * b->inv_mass;
      b->alpha -= r->r.jb[2] * lam * b->inv_inertia;
    }

    if (r->idx == 0) {
      r->el->last_force = lam;
    }
  }
}
