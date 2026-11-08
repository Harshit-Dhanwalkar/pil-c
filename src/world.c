/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* world.c - bodies, shapes, world management, RK4 integration, annotations */
#include "pil_internal.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GROW(arr, n, cap)                                                      \
  do {                                                                         \
    if ((n) >= (cap)) {                                                        \
      (cap) = (cap) ? (cap) * 2 : 8;                                           \
      (arr) = realloc((arr), (size_t)(cap) * sizeof(*(arr)));                  \
      if (!(arr))                                                              \
        abort();                                                               \
    }                                                                          \
  } while (0)

#define TWO_PI 6.28318530717958647692

// shapes
pil_shape pil_point_mass(double r) {
  pil_shape s;
  memset(&s, 0, sizeof s);
  s.kind = PIL_SHAPE_POINT;
  s.r = r;

  return s;
}

pil_shape pil_rect(double w, double h) {
  pil_shape s;
  memset(&s, 0, sizeof s);
  s.kind = PIL_SHAPE_RECT;
  s.w = w;
  s.h = h;

  return s;
}

pil_shape pil_disc(double r) {
  pil_shape s;
  memset(&s, 0, sizeof s);
  s.kind = PIL_SHAPE_DISC;
  s.r = r;

  return s;
}

pil_shape pil_poly(const pil_vec2 *pts, int n) {
  pil_shape s;
  memset(&s, 0, sizeof s);
  s.kind = PIL_SHAPE_POLY;
  s.pts = pts;
  s.npts = n;

  return s;
}

pil_shape pil_arc_band(double r_in, double r_out, double a0, double a1) {
  pil_shape s;
  memset(&s, 0, sizeof s);
  s.kind = PIL_SHAPE_ARC_BAND;
  s.r_in = r_in;
  s.r_out = r_out;
  s.a0 = a0;
  s.a1 = a1;

  return s;
}

// Inertia about COM for unit-agnostic mass m; may also shift s->offset so that
// body origin coincides with centroid
static double shape_inertia(pil_shape *s, double m, int auto_center) {
  switch (s->kind) {
  case PIL_SHAPE_RECT: {
    return m * (s->w * s->w + s->h * s->h) / 12.0;
  }
  case PIL_SHAPE_DISC: {
    return 0.5 * m * s->r * s->r;
  }
  case PIL_SHAPE_POINT: {
    double r = s->r > 0.01 ? s->r : 0.01;
    return 0.5 * m * r * r;
  }
  case PIL_SHAPE_POLY: {
    double A = 0;
    double cx = 0;
    double cy = 0;
    double J = 0;
    for (int i = 0; i < s->npts; i++) {
      pil_vec2 p = s->pts[i], q = s->pts[(i + 1) % s->npts];
      double cr = v_cross(p, q);
      A += cr;
      cx += (p.x + q.x) * cr;
      cy += (p.y + q.y) * cr;
      J += cr * (p.x * p.x + p.x * q.x + q.x * q.x + p.y * p.y + p.y * q.y +
                 q.y * q.y);
    }
    A *= 0.5;
    J /= 12.0;
    if (fabs(A) < 1e-14) {
      return m * 1e-4;
    }
    cx /= 6.0 * A;
    cy /= 6.0 * A;
    double I = m * J / A - m * (cx * cx + cy * cy);
    if (auto_center) {
      s->offset = pil_v2(-cx, -cy);
    }
    return I;
  }
  case PIL_SHAPE_ARC_BAND: {
    double al = 0.5 * (s->a1 - s->a0);
    double pm = 0.5 * (s->a0 + s->a1);
    double ro = s->r_out;
    double ri = s->r_in;
    double rbar = al > 1e-9
                      ? (2.0 * sin(al) / (3.0 * al)) *
                            (ro * ro * ro - ri * ri * ri) / (ro * ro - ri * ri)
                      : ro;
    double Ic = 0.5 * m * (ro * ro + ri * ri);
    if (auto_center) {
      s->offset = pil_v2(-rbar * cos(pm), -rbar * sin(pm));
    }
    return Ic - m * rbar * rbar;
  }
  }

  return m;
}

// Bodies
pil_body *pil_add_body(pil_world *w, const pil_body_desc *d) {
  pil_body *b = calloc(1, sizeof *b);
  if (!b) {
    abort();
  }

  b->id = w->nbodies;
  b->shape = d->shape;
  if (d->shape.kind == PIL_SHAPE_POLY && d->shape.npts > 0) {
    pil_vec2 *copy = malloc((size_t)d->shape.npts * sizeof *copy);
    if (!copy) {
      abort();
    }
    memcpy(copy, d->shape.pts, (size_t)d->shape.npts * sizeof *copy);
    b->shape.pts = copy;
  }

  b->flags = d->flags;
  b->mass = d->mass > 0 ? d->mass : 1.0;
  if (d->inertia > 0) {
    b->inertia = d->inertia;
  } else {
    b->inertia = shape_inertia(&b->shape, b->mass, 1);
  }
  if (d->inertia > 0) {
    shape_inertia(&b->shape, b->mass, 0);
  }
  if (b->inertia <= 0) {
    b->inertia = 1e-6;
  }
  if (d->flags & PIL_STATIC) {
    b->dof = -1;
    b->inv_mass = 0;
    b->inv_inertia = 0;
  } else {
    b->dof = w->ndyn++;
    b->inv_mass = 1.0 / b->mass;
    b->inv_inertia = 1.0 / b->inertia;
  }

  b->pos = d->pos;
  b->angle = d->angle;
  b->vel = d->vel;
  b->omega = d->omega;
  b->pos0 = b->pos;
  b->angle0 = b->angle;
  b->vel0 = b->vel;
  b->omega0 = b->omega;
  b->fill = d->fill;
  pil_body_set_label(b, d->label, d->label_dx, d->label_dy);
  GROW(w->bodies, w->nbodies, w->bcap);
  w->bodies[w->nbodies++] = b;

  return b;
}

void pil_body_set_label(pil_body *b, const char *text, double dx, double dy) {
  snprintf(b->label, sizeof b->label, "%s", text ? text : "");
  b->label_dx = dx;
  b->label_dy = dy;
}

void pil_body_set_state(pil_body *b, pil_vec2 pos, double angle, pil_vec2 vel,
                        double omega) {
  b->pos = b->pos0 = pos;
  b->angle = b->angle0 = angle;
  b->vel = b->vel0 = vel;
  b->omega = b->omega0 = omega;
}

pil_vec2 pil_body_point(const pil_body *b, pil_vec2 local) {
  return v_add(b->pos, v_rot(local, b->angle));
}

pil_anchor pil_on_at(pil_body *b, double x, double y) {
  pil_vec2 l = v_rot(v_sub(pil_v2(x, y), b->pos), -b->angle);

  return pil_on(b, l.x, l.y);
}

pil_vec2 pil_poly_centroid(const pil_vec2 *pts, int n) {
  double A = 0;
  double cx = 0;
  double cy = 0;
  for (int i = 0; i < n; i++) {
    pil_vec2 p = pts[i];
    pil_vec2 q = pts[(i + 1) % n];
    double cr = v_cross(p, q);
    A += cr;
    cx += (p.x + q.x) * cr;
    cy += (p.y + q.y) * cr;
  }

  if (fabs(A) < 1e-14) {
    return n > 0 ? pts[0] : pil_v2(0, 0);
  }

  return pil_v2(cx / (3.0 * A), cy / (3.0 * A));
}

pil_vec2 pil_body_shape_origin(const pil_body *b) { return b->shape.offset; }

pil_vec2 pil_body_velocity_at(const pil_body *b, pil_vec2 p) {
  return v_add(b->vel, v_mul(v_perp(v_sub(p, b->pos)), b->omega));
}

void pil_apply_force(pil_body *b, pil_vec2 p, pil_vec2 f) {
  if (!b || b->dof < 0) {
    return;
  }

  b->force = v_add(b->force, f);
  b->torque += v_cross(v_sub(p, b->pos), f);
}

void pil_apply_torque(pil_body *b, double tau) {
  if (b && b->dof >= 0) {
    b->torque += tau;
  }
}

// World
pil_world *pil_world_new(void) {
  pil_world *w = calloc(1, sizeof *w);
  if (!w) {
    abort();
  }
  w->g = pil_v2(0, -9.81);
  w->alpha = w->beta = 10.0;

  return w;
}

void pil_world_free(pil_world *w) {
  if (!w) {
    return;
  }

  for (int i = 0; i < w->nbodies; i++) {
    if (w->bodies[i]->shape.kind == PIL_SHAPE_POLY) {
      free((void *)w->bodies[i]->shape.pts);
    }
    free(w->bodies[i]);
  }

  for (int i = 0; i < w->nels; i++) {
    free(w->els[i]);
  }
  for (int i = 0; i < w->ndecos; i++) {
    free(w->decos[i]);
  }
  for (int i = 0; i < w->nterr; i++) {
    free(w->terrains[i]);
  }
  pil_plots_free(w);
  pil_circuits_free(w);
  pil_lsys_free_all(w);
  free(w->bodies);
  free(w->els);
  free(w->decos);
  free(w->terrains);
  free(w->overlays);
  free(w->K);
  free(w->rows);
  free(w->s0);
  free(w->k1);
  free(w->k2);
  free(w->k3);
  free(w->k4);
  free(w->tmp);
  free(w);
}

void pil_world_set_gravity(pil_world *w, double gx, double gy) {
  w->g = pil_v2(gx, gy);
}
void pil_world_set_stabilization(pil_world *w, double a, double b) {
  w->alpha = a;
  w->beta = b;
}
void pil_world_set_drag(pil_world *w, double lin, double ang) {
  w->drag_lin = lin;
  w->drag_ang = ang;
}

double pil_world_time(const pil_world *w) { return w->t; }
int pil_world_body_count(const pil_world *w) { return w->nbodies; }
pil_body *pil_world_body(pil_world *w, int i) {
  return (i >= 0 && i < w->nbodies) ? w->bodies[i] : NULL;
}

void pil_world_add_overlay(pil_world *w, pil_overlay_fn fn, void *user) {
  GROW(w->overlays, w->noverlays, w->ocap);
  w->overlays[w->noverlays].fn = fn;
  w->overlays[w->noverlays].user = user;
  w->noverlays++;
}

// RK4
static void load_state(pil_world *w, const double *s) {
  for (int i = 0; i < w->nbodies; i++) {
    pil_body *b = w->bodies[i];
    if (b->dof < 0) {
      continue;
    }

    const double *p = s + 6 * b->dof;
    b->pos = pil_v2(p[0], p[1]);
    b->angle = p[2];
    b->vel = pil_v2(p[3], p[4]);
    b->omega = p[5];
  }
}

static void deriv(pil_world *w, const double *s, double t, double *ds) {
  load_state(w, s);
  w->tcur = t;
  pil_compute_accel(w);
  for (int i = 0; i < w->nbodies; i++) {
    pil_body *b = w->bodies[i];
    if (b->dof < 0) {
      continue;
    }

    double *d = ds + 6 * b->dof;
    d[0] = b->vel.x;
    d[1] = b->vel.y;
    d[2] = b->omega;
    d[3] = b->acc.x;
    d[4] = b->acc.y;
    d[5] = b->alpha;
  }
}

void pil_world_step(pil_world *w, double dt) {
  int n = 6 * w->ndyn;
  if (n == 0) {
    if (!w->frozen) {
      pil_lsys_step_all(w, w->t, dt);
      w->t += dt;
      pil_circuits_step(w, dt);
    }

    return;
  }

  if (n > w->scap) {
    w->scap = n;
    w->s0 = realloc(w->s0, n * sizeof(double));
    w->k1 = realloc(w->k1, n * sizeof(double));
    w->k2 = realloc(w->k2, n * sizeof(double));
    w->k3 = realloc(w->k3, n * sizeof(double));
    w->k4 = realloc(w->k4, n * sizeof(double));
    w->tmp = realloc(w->tmp, n * sizeof(double));
    if (!w->s0 || !w->k1 || !w->k2 || !w->k3 || !w->k4 || !w->tmp) {
      abort();
    }
  }

  for (int i = 0; i < w->nbodies; i++) {
    pil_body *b = w->bodies[i];
    if (b->dof < 0) {
      continue;
    }

    double *p = w->s0 + 6 * b->dof;
    p[0] = b->pos.x;
    p[1] = b->pos.y;
    p[2] = b->angle;
    p[3] = b->vel.x;
    p[4] = b->vel.y;
    p[5] = b->omega;
  }

  double t = w->t;
  double h = w->frozen ? 0.0 : 1.0;
  deriv(w, w->s0, t, w->k1);
  for (int i = 0; i < n; i++) {
    w->tmp[i] = w->s0[i] + 0.5 * dt * w->k1[i];
  }

  deriv(w, w->tmp, t + 0.5 * dt * h, w->k2);
  for (int i = 0; i < n; i++) {
    w->tmp[i] = w->s0[i] + 0.5 * dt * w->k2[i];
  }

  deriv(w, w->tmp, t + 0.5 * dt * h, w->k3);
  for (int i = 0; i < n; i++) {
    w->tmp[i] = w->s0[i] + dt * w->k3[i];
  }

  deriv(w, w->tmp, t + dt * h, w->k4);
  for (int i = 0; i < n; i++) {
    w->tmp[i] = w->s0[i] +
                dt / 6.0 * (w->k1[i] + 2 * w->k2[i] + 2 * w->k3[i] + w->k4[i]);
  }

  load_state(w, w->tmp);
  if (!w->frozen) {
    pil_lsys_step_all(w, w->t, dt);
    w->t += dt;
    pil_circuits_step(w, dt);
  }

  w->tcur = w->t;

  pil_project_rope_snaps(w);
}

void pil_world_advance(pil_world *w, double duration, double max_dt) {
  if (duration <= 0) {
    return;
  }

  int n = (int)ceil(duration / max_dt);
  if (n < 1) {
    n = 1;
  }

  double dt = duration / n;
  for (int i = 0; i < n; i++) {
    pil_world_step(w, dt);
  }
}

// Diagnostics
double pil_world_kinetic(const pil_world *w) {
  double T = 0;
  for (int i = 0; i < w->nbodies; i++) {
    const pil_body *b = w->bodies[i];
    if (b->dof >= 0) {
      T += 0.5 * b->mass * v_dot(b->vel, b->vel) +
           0.5 * b->inertia * b->omega * b->omega;
    }
  }

  return T;
}

double pil_world_dissipated(const pil_world *w) { return w->dissipated; }

double pil_world_energy(const pil_world *w) {
  double E = 0;
  for (int i = 0; i < w->nbodies; i++) {
    const pil_body *b = w->bodies[i];
    if (b->dof < 0) {
      continue;
    }
    E += 0.5 * b->mass * v_dot(b->vel, b->vel) +
         0.5 * b->inertia * b->omega * b->omega;
    E -= b->mass * v_dot(w->g, b->pos);
  }

  for (int i = 0; i < w->nels; i++) {
    E += pil_element_energy(w->els[i], w->t);
  }

  return E;
}

double pil_world_constraint_error(const pil_world *w) {
  double m = 0;
  for (int i = 0; i < w->nels; i++) {
    const pil_element *e = w->els[i];
    if (!pil_is_constraint(e)) {
      continue;
    }

    pil_row r[2];
    int n = pil_constraint_eval(e, w->t, r, 1);
    for (int k = 0; k < n; k++) {
      double v = r[k].phi;
      if (pil_is_rope(e) && v < 0) {
        v = 0; // slack is fine
      }
      if (fabs(v) > m) {
        m = fabs(v);
      }
    }
  }
  return m;
}

// snapshot / reset / settle
void pil_world_snapshot(pil_world *w) {
  for (int i = 0; i < w->nbodies; i++) {
    pil_body *b = w->bodies[i];
    b->pos0 = b->pos;
    b->angle0 = b->angle;
    b->vel0 = b->vel;
    b->omega0 = b->omega;
  }
  w->have_snapshot = 1;
}

void pil_world_reset(pil_world *w) {
  for (int i = 0; i < w->nbodies; i++) {
    pil_body *b = w->bodies[i];
    b->pos = b->pos0;
    b->angle = b->angle0;
    b->vel = b->vel0;
    b->omega = b->omega0;
  }

  w->t = 0;
  w->tcur = 0;
  w->dissipated = 0;
  for (int i = 0; i < w->nels; i++) {
    w->els[i]->rope_off = 0;
  }

  pil_circuits_reset(w);
  pil_elements_rebase(w);
  pil_lsys_reset_all(w);
}

void pil_world_settle(pil_world *w, double seconds) {
  double dl = w->drag_lin;
  double da = w->drag_ang;
  int fr = w->frozen;
  w->drag_lin = 4.0;
  w->drag_ang = 4.0;
  w->frozen = 1;
  int n = (int)(seconds / 2e-4);
  for (int i = 0; i < n; i++) {
    pil_world_step(w, 2e-4);
  }

  w->drag_lin = dl;
  w->drag_ang = da;
  w->frozen = fr;
  for (int i = 0; i < w->nbodies; i++) {
    w->bodies[i]->vel = pil_v2(0, 0);
    w->bodies[i]->omega = 0;
  }

  for (int i = 0; i < w->nels; i++) {
    if (w->els[i]->kind == EL_ROT_DRIVE) {
      w->els[i]->a.body->omega = w->els[i]->omega;
    }
  }

  w->t = 0;

  pil_world_snapshot(w);
  pil_elements_rebase(w);
}

// picking / grabbing
static int in_shape(const pil_shape *s, pil_vec2 q) {
  switch (s->kind) {
  case PIL_SHAPE_RECT: {
    return fabs(q.x) <= 0.5 * s->w && fabs(q.y) <= 0.5 * s->h;
  }
  case PIL_SHAPE_DISC:
  case PIL_SHAPE_POINT: {
    double r = s->r > 0.03 ? s->r : 0.03;
    return v_dot(q, q) <= r * r;
  }
  case PIL_SHAPE_POLY: {
    int in = 0;
    for (int i = 0, j = s->npts - 1; i < s->npts; j = i++) {
      pil_vec2 a = s->pts[i], b = s->pts[j];
      if ((a.y > q.y) != (b.y > q.y) &&
          q.x < (b.x - a.x) * (q.y - a.y) / (b.y - a.y) + a.x) {
        in = !in;
      }
    }

    return in;
  }
  case PIL_SHAPE_ARC_BAND: {
    double r = sqrt(v_dot(q, q));
    if (r < s->r_in || r > s->r_out) {
      return 0;
    }
    double da = fmod(atan2(q.y, q.x) - s->a0, TWO_PI);
    if (da < 0) {
      da += TWO_PI;
    }
    return da <= s->a1 - s->a0;
  }
  }

  return 0;
}

pil_body *pil_world_pick(const pil_world *w, pil_vec2 pt, pil_vec2 *local_out) {
  for (int i = w->nbodies - 1; i >= 0; i--) {
    pil_body *b = w->bodies[i];
    if (b->dof < 0 || (b->flags & PIL_HIDDEN)) {
      continue;
    }

    pil_vec2 loc = v_rot(v_sub(pt, b->pos), -b->angle);
    if (in_shape(&b->shape, v_sub(loc, b->shape.offset))) {
      if (local_out) {
        *local_out = loc;
      }
      return b;
    }
  }

  return NULL;
}

void pil_world_set_grab(pil_world *w, pil_body *b, pil_vec2 local,
                        pil_vec2 target, double k, double c) {
  w->grab.b = b;
  w->grab.local = local;
  w->grab.target = target;
  w->grab.k = k;
  w->grab.c = c;
}

void pil_world_clear_grab(pil_world *w) { w->grab.b = NULL; }

// annotations
static pil_deco *new_deco(pil_world *w, pil_deco_kind k, const char *text) {
  pil_deco *d = calloc(1, sizeof *d);
  if (!d) {
    abort();
  }

  d->kind = k;
  snprintf(d->text, sizeof d->text, "%s", text ? text : "");
  GROW(w->decos, w->ndecos, w->dcap);
  w->decos[w->ndecos++] = d;

  return d;
}

pil_deco *pil_add_energy_bar(pil_world *w, pil_vec2 bl, double w_px,
                             double h_px, double v_ref) {
  pil_deco *d = new_deco(w, D_ENERGY_BAR, "");
  d->a = pil_at(bl.x, bl.y);
  d->dx = w_px;
  d->dy = h_px;
  d->ref = v_ref;
  d->offset = pil_world_energy(w) - v_ref; // E0 above potential zero

  return d;
}

pil_deco *pil_add_gravity_into_page(pil_world *w, pil_vec2 at) {
  pil_deco *d = pil_add_gravity_indicator(w, at);
  d->ref = 1;

  return d;
}

pil_deco *pil_add_label(pil_world *w, pil_anchor at, const char *text,
                        double dx, double dy) {
  pil_deco *d = new_deco(w, D_LABEL, text);
  d->a = at;
  d->dx = dx;
  d->dy = dy;

  return d;
}

pil_deco *pil_add_dimension(pil_world *w, pil_anchor a, pil_anchor b,
                            double offset, const char *text) {
  pil_deco *d = new_deco(w, D_DIM, text);
  d->a = a;
  d->b = b;
  d->offset = offset;

  return d;
}

pil_deco *pil_add_leader(pil_world *w, pil_anchor at, double dx, double dy,
                         const char *text) {
  pil_deco *d = new_deco(w, D_ARROW, text);
  d->a = at;
  d->dx = dx;
  d->dy = dy;
  d->reverse = 1;

  return d;
}

pil_deco *pil_add_angle_between(pil_world *w, pil_anchor c, pil_anchor p1,
                                pil_anchor p2, double radius_px,
                                const char *text) {
  pil_deco *d = new_deco(w, D_ANGLE_PTS, text);
  d->a = c;
  d->b = p1;
  d->c = p2;
  d->radius = radius_px;

  return d;
}

pil_deco *pil_add_arrow(pil_world *w, pil_anchor from, double dx, double dy,
                        const char *text) {
  pil_deco *d = new_deco(w, D_ARROW, text);
  d->a = from;
  d->dx = dx;
  d->dy = dy;

  return d;
}

pil_deco *pil_add_line(pil_world *w, pil_anchor a, pil_anchor b, int dash_dot) {
  pil_deco *d = new_deco(w, D_LINE, "");
  d->a = a;
  d->b = b;
  d->dashed = dash_dot;

  return d;
}

pil_deco *pil_add_ground(pil_world *w, pil_vec2 a, pil_vec2 b) {
  pil_deco *d = new_deco(w, D_GROUND, "");
  d->a = pil_at(a.x, a.y);
  d->b = pil_at(b.x, b.y);

  return d;
}

pil_deco *pil_add_gravity_indicator(pil_world *w, pil_vec2 at) {
  pil_deco *d = new_deco(w, D_GRAVITY, "g");
  d->a = pil_at(at.x, at.y);

  return d;
}

pil_deco *pil_add_angle_marker(pil_world *w, pil_anchor c, double ref,
                               pil_body *body, double radius_px,
                               const char *text) {
  pil_deco *d = new_deco(w, D_ANGLE, text);
  d->a = c;
  d->ref = ref;
  d->body = body;
  d->radius = radius_px;

  return d;
}

void pil_deco_hide(pil_deco *d) { d->hidden = 1; }
