/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* elements.c - force elements and constraints (things to "drop in") */

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

// Anchors
pil_vec2 pil_anchor_world(pil_anchor a) {
  return a.body ? v_add(a.body->pos, v_rot(a.local, a.body->angle)) : a.local;
}

pil_vec2 pil_anchor_r(pil_anchor a) {
  return a.body ? v_rot(a.local, a.body->angle) : pil_v2(0, 0);
}

pil_vec2 pil_anchor_vel(pil_anchor a) {
  if (!a.body) {
    return pil_v2(0, 0);
  }

  return v_add(a.body->vel,
               v_mul(v_perp(v_rot(a.local, a.body->angle)), a.body->omega));
}

static double a_omega(pil_anchor a) { return a.body ? a.body->omega : 0.0; }

// one coordinate (x, y, angle) of a body and its rate; world has 0
static double coord_q(const pil_body *b, int c) {
  if (!b) {
    return 0.0;
  }

  return c == PIL_COORD_X ? b->pos.x : (c == PIL_COORD_Y ? b->pos.y : b->angle);
}

static double coord_v(const pil_body *b, int c) {
  if (!b) {
    return 0.0;
  }

  return c == PIL_COORD_X ? b->vel.x : (c == PIL_COORD_Y ? b->vel.y : b->omega);
}

// Creation
static pil_element *new_el(pil_world *w, pil_el_kind k) {
  pil_element *e = calloc(1, sizeof *e);
  if (!e) {
    abort();
  }

  e->world = w;
  e->kind = k;
  e->label_dx = -28;
  e->label_dy = 0;
  GROW(w->els, w->nels, w->ecap);
  w->els[w->nels++] = e;

  return e;
}

pil_element *pil_add_link(pil_world *w, pil_body *a, pil_coord ca, pil_body *b,
                          pil_coord cb, double ratio) {
  pil_element *e = new_el(w, EL_LINK);
  e->a.body = a;
  e->b.body = b;
  e->ca = ca;
  e->cb = cb;
  e->k = ratio;
  e->phase = coord_q(a, ca) - ratio * coord_q(b, cb);
  e->hidden = 1;

  return e;
}

static double anchor_dist(pil_anchor a, pil_anchor b) {
  pil_vec2 d = v_sub(pil_anchor_world(a), pil_anchor_world(b));

  return sqrt(v_dot(d, d));
}

pil_element *pil_add_spring(pil_world *w, pil_anchor a, pil_anchor b, double k,
                            double rest) {
  pil_element *e = new_el(w, EL_SPRING);
  e->a = a;
  e->b = b;
  e->k = k;
  e->length = anchor_dist(a, b);
  e->rest = rest < 0 ? e->length : rest;

  return e;
}

pil_element *pil_add_damper(pil_world *w, pil_anchor a, pil_anchor b,
                            double c) {
  pil_element *e = new_el(w, EL_DAMPER);
  e->a = a;
  e->b = b;
  e->c = c;
  e->length = anchor_dist(a, b);

  return e;
}

pil_element *pil_add_torsion_spring(pil_world *w, pil_body *a, pil_body *b,
                                    double k, double rest_angle, double c) {
  pil_element *e = new_el(w, EL_TORSION);
  e->a.body = a;
  e->b.body = b;
  e->k = k;
  e->rest = rest_angle;
  e->c = c;

  return e;
}

pil_element *pil_add_force(pil_world *w, pil_anchor at, pil_vec2 f) {
  pil_element *e = new_el(w, EL_FORCE);
  e->a = at;
  e->f = f;

  return e;
}

pil_element *pil_add_custom_force(pil_world *w, pil_force_fn fn, void *user) {
  pil_element *e = new_el(w, EL_CUSTOM);
  e->fn = fn;
  e->user = user;
  e->hidden = 1;

  return e;
}

pil_element *pil_add_revolute(pil_world *w, pil_anchor a, pil_anchor b) {
  pil_element *e = new_el(w, EL_REVOLUTE);
  e->a = a;
  e->b = b;

  return e;
}

pil_element *pil_add_distance(pil_world *w, pil_anchor a, pil_anchor b,
                              double length, unsigned flags) {
  pil_element *e = new_el(w, EL_DISTANCE);
  e->a = a;
  e->b = b;
  e->flags = flags;
  e->length = length < 0 ? anchor_dist(a, b) : length;

  return e;
}

static pil_vec2 unit(pil_vec2 d) {
  double l = sqrt(v_dot(d, d));

  return l > 0 ? v_mul(d, 1.0 / l) : pil_v2(1, 0);
}

int pil_is_rope(const pil_element *e) {
  return e->kind == EL_PULLEY ||
         (e->kind == EL_DISTANCE && (e->flags & PIL_ROPE));
}

double pil_rope_length(const pil_element *e) {
  pil_vec2 pa = pil_anchor_world(e->a), pb = pil_anchor_world(e->b);
  if (e->kind == EL_PULLEY) {
    pil_vec2 da = v_sub(pa, e->p1), db = v_sub(pb, e->p2);

    return sqrt(v_dot(da, da)) + sqrt(v_dot(db, db));
  }

  pil_vec2 d = v_sub(pa, pb);

  return sqrt(v_dot(d, d));
}

pil_element *pil_add_pulley_rope(pil_world *w, pil_anchor a, pil_anchor b,
                                 pil_vec2 c1, pil_vec2 c2, double length) {
  pil_element *e = new_el(w, EL_PULLEY);
  e->a = a;
  e->b = b;
  e->p1 = c1;
  e->p2 = c2;
  e->length = length < 0 ? pil_rope_length(e) : length;

  return e;
}

pil_element *pil_add_slot(pil_world *w, pil_anchor a, pil_anchor line_origin,
                          pil_vec2 dir) {
  pil_element *e = new_el(w, EL_SLOT);
  e->a = a;
  e->b = line_origin;
  e->dir = unit(dir);

  return e;
}

pil_element *pil_add_fixed_angle(pil_world *w, pil_body *a, pil_body *b,
                                 double rel_angle) {
  pil_element *e = new_el(w, EL_FIXED_ANGLE);
  e->a.body = a;
  e->b.body = b;
  e->phase = rel_angle;

  return e;
}

pil_element *pil_add_stop(pil_world *w, pil_body *b, pil_coord coord,
                          double limit, int upper, double restitution) {
  pil_element *e = new_el(w, EL_STOP);
  e->a.body = b;
  e->ca = coord;
  e->phase = limit;
  e->sigma = upper ? 1 : -1;
  e->k = restitution;
  e->hidden = 1;

  return e;
}

double pil_stop_row(const pil_element *e, pil_row *r) {
  memset(r, 0, sizeof *r);
  r->phi = e->sigma * (coord_q(e->a.body, e->ca) - e->phase);
  r->vel = e->sigma * coord_v(e->a.body, e->ca);
  r->ja[e->ca] = e->sigma;

  return e->k;
}

pil_element *pil_add_prismatic(pil_world *w, pil_anchor a,
                               pil_anchor line_origin, pil_vec2 dir) {
  double ta = a.body ? a.body->angle : 0;
  double tb = line_origin.body ? line_origin.body->angle : 0;
  pil_element *s = pil_add_slot(w, a, line_origin, dir);
  if (a.body) {
    pil_add_fixed_angle(w, a.body, line_origin.body, ta - tb)->hidden = 1;
  }

  return s;
}

pil_element *pil_add_rotation_drive(pil_world *w, pil_body *a, double omega) {
  pil_element *e = new_el(w, EL_ROT_DRIVE);
  e->a.body = a;
  e->omega = omega;
  e->phase = a->angle;
  e->hidden = 1;

  return e;
}

pil_element *pil_add_rolling(pil_world *w, pil_body *b, pil_vec2 center_local,
                             double radius, pil_vec2 line_point,
                             pil_vec2 line_dir) {
  pil_element *e = new_el(w, EL_ROLLING);
  e->a.body = b;
  e->a.local = center_local;
  e->radius = radius;
  e->b = pil_at(line_point.x, line_point.y);
  e->dir = unit(line_dir);
  e->hidden = 1;

  pil_elements_rebase(w);

  return e;
}

pil_element *pil_add_circle_rolling(pil_world *w, pil_anchor a, double ra,
                                    pil_anchor b, double rb, int inside) {
  pil_element *e = new_el(w, EL_CIRCLE_ROLL);
  e->a = a;
  e->b = b;
  e->radius = ra;
  e->radius2 = rb;
  e->sigma = inside ? -1 : +1;
  e->hidden = 1;

  pil_elements_rebase(w);

  return e;
}

pil_element *pil_add_weld(pil_world *w, pil_anchor a, pil_anchor b) {
  pil_element *pin = pil_add_revolute(w, a, b);
  double ta = a.body ? a.body->angle : 0, tb = b.body ? b.body->angle : 0;
  if (a.body) {
    pil_add_fixed_angle(w, a.body, b.body, ta - tb)->hidden = 1;
  } else if (b.body) {
    pil_add_fixed_angle(w, b.body, a.body, tb - ta)->hidden = 1;
  }
  pin->hidden = 1;

  return pin;
}

pil_element *pil_add_floor_contact(pil_world *w, pil_body *body,
                                   pil_anchor line_point, pil_vec2 line_dir,
                                   double k, double c, double mu) {
  pil_element *e = new_el(w, EL_FLOOR);
  e->a.body = body;
  e->b = line_point;
  e->dir = unit(line_dir);
  e->k = k;
  e->c = c;
  e->mu = mu;
  e->hidden = 1;

  return e;
}

void pil_floor_set_range(pil_element *e, double s0, double s1) {
  e->ranged = 1;
  e->s0 = s0;
  e->s1 = s1;
}

pil_element *pil_add_circle_contact(pil_world *w, pil_body *body,
                                    pil_anchor centre, double R, int inside,
                                    double k, double c, double mu) {
  pil_element *e = new_el(w, EL_CIRCLE_CONTACT);
  e->a.body = body;
  e->b = centre;
  e->radius2 = R;
  e->sigma = inside ? -1 : 1;
  e->k = k;
  e->c = c;
  e->mu = mu;
  e->hidden = 1;

  return e;
}

void pil_circle_contact_set_arc(pil_element *e, double a0, double span) {
  e->ranged = 1;
  e->s0 = a0;
  e->s1 = span;
}

pil_terrain *pil_add_terrain(pil_world *w, pil_height_fn fn, void *user,
                             double scroll_speed) {
  pil_terrain *t = calloc(1, sizeof *t);
  if (!t) {
    abort();
  }

  t->fn = fn;
  t->user = user;
  t->scroll = scroll_speed;
  GROW(w->terrains, w->nterr, w->tcap);
  w->terrains[w->nterr++] = t;

  return t;
}

double pil_terrain_height(const pil_terrain *t, double x, double time) {
  return t->fn(x + t->scroll * time, t->user);
}

pil_element *pil_add_terrain_contact(pil_world *w, pil_terrain *t,
                                     pil_body *wheel, double radius, double k,
                                     double c) {
  pil_element *e = new_el(w, EL_TERRAIN_CONTACT);
  e->terrain = t;
  e->a.body = wheel;
  e->radius = radius;
  e->k = k;
  e->c = c;
  e->hidden = 1;

  return e;
}

// Generic property setters
void pil_element_set_label(pil_element *e, const char *text) {
  snprintf(e->label, sizeof e->label, "%s", text ? text : "");
}

void pil_element_set_label_at(pil_element *e, const char *text, double dx,
                              double dy) {
  pil_element_set_label(e, text);
  e->label_dx = dx;
  e->label_dy = dy;
}

void pil_element_set_drawing(pil_element *e, int coils, double width_px) {
  e->coils = coils;
  e->width_px = width_px;
}

void pil_element_hide(pil_element *e) { e->hidden = 1; }
double pil_element_force(const pil_element *e) { return e->last_force; }

// Rebase: phases that depend on current configuration
void pil_elements_rebase(pil_world *w) {
  for (int i = 0; i < w->nels; i++) {
    pil_element *e = w->els[i];
    if (e->kind == EL_ROT_DRIVE) {
      e->phase = e->a.body->angle - e->omega * w->t;
    } else if (e->kind == EL_CIRCLE_ROLL) {
      pil_vec2 d = v_sub(pil_anchor_world(e->a), pil_anchor_world(e->b));
      double thb = e->b.body ? e->b.body->angle : 0;
      double thr = e->a.body ? e->a.body->angle : 0;
      pil_vec2 dl = v_rot(d, -thb);
      double D = e->radius2 + e->sigma * e->radius;
      e->ref = atan2(dl.y, dl.x);
      e->phase = D * e->ref + e->sigma * e->radius * (thb - thr);
    } else if (e->kind == EL_LINK) {
      e->phase = coord_q(e->a.body, e->ca) - e->k * coord_q(e->b.body, e->cb);
    } else if (e->kind == EL_ROLLING) {
      pil_body *b = e->a.body;
      pil_vec2 C = pil_anchor_world(e->a);
      e->phase = v_dot(e->dir, v_sub(C, e->b.local)) + e->radius * b->angle;
    }
  }
}

// Floor contact helper
static int contact_points(const pil_body *b, pil_vec2 n, pil_vec2 *pts) {
  const pil_shape *s = &b->shape;
  pil_vec2 c0 = v_add(b->pos, v_rot(s->offset, b->angle));
  int np = 0;
  switch (s->kind) {
  case PIL_SHAPE_POINT:
  case PIL_SHAPE_DISC:
    pts[np++] = v_sub(c0, v_mul(n, s->r));
    break;
  case PIL_SHAPE_RECT: {
    for (int i = 0; i < 4; i++) {
      pil_vec2 q =
          pil_v2((i & 1 ? 0.5 : -0.5) * s->w, (i & 2 ? 0.5 : -0.5) * s->h);
      pts[np++] = v_add(c0, v_rot(q, b->angle));
    }
    break;
  }
  case PIL_SHAPE_POLY:
    for (int i = 0; i < s->npts && np < 64; i++) {
      pts[np++] = v_add(c0, v_rot(s->pts[i], b->angle));
    }
    break;
  default:
    break;
  }

  return np;
}

static double floor_eval(const pil_element *e, int apply) {
  pil_body *bd = e->a.body;
  pil_anchor L = e->b;
  pil_vec2 P0 = pil_anchor_world(L);
  pil_vec2 u = v_rot(e->dir, L.body ? L.body->angle : 0.0);
  pil_vec2 n = v_perp(u);
  pil_vec2 pts[64];

  int np = contact_points(bd, n, pts);
  double energy = 0;
  double vs = 1e-2;
  for (int i = 0; i < np; i++) {
    double depth = -v_dot(n, v_sub(pts[i], P0));
    if (depth <= 0) {
      continue;
    }

    if (e->ranged) {
      double s = v_dot(u, v_sub(pts[i], P0));
      if (s < e->s0 || s > e->s1) {
        continue;
      }
    }

    energy += 0.5 * e->k * depth * depth;
    if (!apply) {
      continue;
    }

    pil_vec2 vp = pil_body_velocity_at(bd, pts[i]);
    if (L.body) {
      vp = v_sub(vp, pil_body_velocity_at(L.body, pts[i]));
    }

    double N = e->k * depth - e->c * v_dot(n, vp);
    if (N < 0) {
      N = 0;
    }

    double ft = -v_dot(u, vp) * (e->mu * N / vs), lim = e->mu * N;
    if (ft > lim) {
      ft = lim;
    } else if (ft < -lim) {
      ft = -lim;
    }

    pil_vec2 F = v_add(v_mul(n, N), v_mul(u, ft));
    pil_apply_force(bd, pts[i], F);
    if (L.body) {
      pil_apply_force(L.body, pts[i], v_mul(F, -1.0));
    }
  }

  return energy;
}

// Disc on a circular track: depth is how far disc has pushed into wall, n
// points from wall towards disc
static double circle_eval(const pil_element *e, int apply) {
  pil_body *bd = e->a.body;
  pil_anchor L = e->b;
  double R = e->radius2, r = bd->shape.r;
  pil_vec2 O = pil_anchor_world(L);
  pil_vec2 C = v_add(bd->pos, v_rot(bd->shape.offset, bd->angle));
  pil_vec2 d = v_sub(C, O);
  double len = sqrt(v_dot(d, d));
  if (len < 1e-9) {
    return 0;
  }

  pil_vec2 out = v_mul(d, 1.0 / len); // from track centre to disc
  double depth = e->sigma < 0 ? len + r - R : R + r - len;
  // disc centre must be on the proper side of the wall: a disc outside a loop
  // is not pushed into it
  if (depth <= 0 || (e->sigma < 0 ? len > R : len < R)) {
    return 0;
  }

  if (e->ranged) {
    double span = e->s1, a = atan2(d.y, d.x) - e->s0;
    a -= 2 * M_PI * floor(a / (2 * M_PI));
    if (a > span) {
      return 0;
    }
  }

  double energy = 0.5 * e->k * depth * depth;
  if (!apply) {
    return energy;
  }

  pil_vec2 n = e->sigma < 0 ? v_mul(out, -1.0) : out;
  pil_vec2 u = v_perp(n), pt = v_sub(C, v_mul(n, r)), vs_ = {0, 0};
  pil_vec2 vp = pil_body_velocity_at(bd, pt);
  if (L.body) {
    vs_ = pil_body_velocity_at(L.body, pt);
    vp = v_sub(vp, vs_);
  }

  double N = e->k * depth - e->c * v_dot(n, vp);
  if (N < 0) {
    N = 0;
  }

  double vsl = 1e-2, ft = -v_dot(u, vp) * (e->mu * N / vsl), lim = e->mu * N;
  if (ft > lim) {
    ft = lim;
  } else if (ft < -lim) {
    ft = -lim;
  }

  pil_vec2 F = v_add(v_mul(n, N), v_mul(u, ft));
  pil_apply_force(bd, pt, F);
  if (L.body) {
    pil_apply_force(L.body, pt, v_mul(F, -1.0));
  }

  return energy;
}

// Force elements
void pil_element_apply(pil_element *e, pil_world *w, double t) {
  switch (e->kind) {
  case EL_SPRING:
  case EL_DAMPER: {
    pil_vec2 pa = pil_anchor_world(e->a), pb = pil_anchor_world(e->b);
    pil_vec2 d = v_sub(pa, pb);
    double len = sqrt(v_dot(d, d));
    if (len < 1e-12) {
      e->last_force = 0;
      break;
    }

    pil_vec2 n = v_mul(d, 1.0 / len);
    double mag; // tension-positive scalar
    if (e->kind == EL_SPRING) {
      mag = e->k * (len - e->rest);
    } else {
      mag = e->c * v_dot(n, v_sub(pil_anchor_vel(e->a), pil_anchor_vel(e->b)));
    }

    e->last_force = mag;
    pil_apply_force(e->a.body, pa, v_mul(n, -mag));
    pil_apply_force(e->b.body, pb, v_mul(n, mag));

    break;
  }
  case EL_TORSION: {
    double ta = e->a.body ? e->a.body->angle : 0;
    double tb = e->b.body ? e->b.body->angle : 0;
    double wa = a_omega(e->a), wb = a_omega(e->b);
    double tau = -e->k * (ta - tb - e->rest) - e->c * (wa - wb);
    e->last_force = tau;
    pil_apply_torque(e->a.body, tau);
    pil_apply_torque(e->b.body, -tau);

    break;
  }
  case EL_FORCE: {
    pil_apply_force(e->a.body, pil_anchor_world(e->a), e->f);

    break;
  }
  case EL_CUSTOM: {
    e->fn(w, t, e->user);

    break;
  }
  case EL_TERRAIN_CONTACT: {
    pil_body *b = e->a.body;
    double ds = 1e-4, s = b->pos.x + e->terrain->scroll * t;
    double h = e->terrain->fn(s, e->terrain->user);
    double hd = e->terrain->scroll *
                (e->terrain->fn(s + ds, e->terrain->user) -
                 e->terrain->fn(s - ds, e->terrain->user)) /
                (2 * ds);
    double delta = h + e->radius - b->pos.y;
    e->last_force = 0;
    if (delta > 0) {
      double F = e->k * delta + e->c * (hd - b->vel.y);
      if (F < 0) {
        F = 0;
      }
      e->last_force = F;
      pil_apply_force(b, b->pos, pil_v2(0, F));
    }

    break;
  }
  case EL_FLOOR: {
    floor_eval(e, 1);

    break;
  }
  case EL_CIRCLE_CONTACT: {
    circle_eval(e, 1);

    break;
  }
  default: {

    break;
  }
  }
}

double pil_element_energy(const pil_element *e, double t) {
  switch (e->kind) {
  case EL_SPRING: {
    pil_vec2 d = v_sub(pil_anchor_world(e->a), pil_anchor_world(e->b));
    double x = sqrt(v_dot(d, d)) - e->rest;

    return 0.5 * e->k * x * x;
  }
  case EL_TORSION: {
    double ta = e->a.body ? e->a.body->angle : 0;
    double tb = e->b.body ? e->b.body->angle : 0;
    double x = ta - tb - e->rest;

    return 0.5 * e->k * x * x;
  }
  case EL_TERRAIN_CONTACT: {
    pil_body *b = e->a.body;
    double h =
        e->terrain->fn(b->pos.x + e->terrain->scroll * t, e->terrain->user);
    double delta = h + e->radius - b->pos.y;

    return delta > 0 ? 0.5 * e->k * delta * delta : 0.0;
  }
  case EL_FLOOR: {
    return floor_eval(e, 0);
  }
  case EL_CIRCLE_CONTACT: {
    return circle_eval(e, 0);
  }
  default:
    return 0.0;
  }
}

// Constraints
int pil_is_constraint(const pil_element *e) {
  return e->kind >= EL_REVOLUTE && e->kind <= EL_PULLEY;
}

static void finish(const pil_element *e, pil_row *r, double g0) {
  r->gamma = g0 - 2.0 * e->world->alpha * r->vel -
             e->world->beta * e->world->beta * r->phi;
}

int pil_constraint_eval(const pil_element *e, double t, pil_row rows[2],
                        int force_active) {
  memset(rows, 0, 2 * sizeof(pil_row));
  switch (e->kind) {

  case EL_REVOLUTE: {
    pil_vec2 ra = pil_anchor_r(e->a), rb = pil_anchor_r(e->b);
    pil_vec2 pa = pil_anchor_world(e->a), pb = pil_anchor_world(e->b);
    pil_vec2 va = pil_anchor_vel(e->a), vb = pil_anchor_vel(e->b);
    double wa2 = a_omega(e->a) * a_omega(e->a);
    double wb2 = a_omega(e->b) * a_omega(e->b);
    rows[0].phi = pa.x - pb.x;
    rows[0].vel = va.x - vb.x;
    rows[0].ja[0] = 1;
    rows[0].ja[2] = -ra.y;
    rows[0].jb[0] = -1;
    rows[0].jb[2] = rb.y;
    finish(e, &rows[0], wa2 * ra.x - wb2 * rb.x);
    rows[1].phi = pa.y - pb.y;
    rows[1].vel = va.y - vb.y;
    rows[1].ja[1] = 1;
    rows[1].ja[2] = ra.x;
    rows[1].jb[1] = -1;
    rows[1].jb[2] = -rb.x;
    finish(e, &rows[1], wa2 * ra.y - wb2 * rb.y);

    return 2;
  }

  case EL_DISTANCE: {
    if ((e->flags & PIL_ROPE) && e->rope_off && !force_active) {
      return 0;
    }

    pil_vec2 ra = pil_anchor_r(e->a), rb = pil_anchor_r(e->b);
    pil_vec2 d = v_sub(pil_anchor_world(e->a), pil_anchor_world(e->b));
    pil_vec2 dd = v_sub(pil_anchor_vel(e->a), pil_anchor_vel(e->b));
    double len = sqrt(v_dot(d, d));
    pil_vec2 n = len > 1e-9 ? v_mul(d, 1.0 / len) : pil_v2(1, 0);
    double wa2 = a_omega(e->a) * a_omega(e->a);
    double wb2 = a_omega(e->b) * a_omega(e->b);
    double ndd = v_dot(n, dd);
    double curv = len > 1e-9 ? (v_dot(dd, dd) - ndd * ndd) / len : 0.0;
    rows[0].phi = len - e->length;
    rows[0].vel = ndd;
    rows[0].ja[0] = n.x;
    rows[0].ja[1] = n.y;
    rows[0].ja[2] = v_cross(ra, n);
    rows[0].jb[0] = -n.x;
    rows[0].jb[1] = -n.y;
    rows[0].jb[2] = -v_cross(rb, n);
    finish(e, &rows[0], v_dot(n, v_sub(v_mul(ra, wa2), v_mul(rb, wb2))) - curv);

    return 1;
  }

  case EL_CIRCLE_ROLL: {
    pil_vec2 rr = pil_anchor_r(e->a);
    pil_vec2 rb = pil_anchor_r(e->b);
    pil_vec2 d = v_sub(pil_anchor_world(e->a), pil_anchor_world(e->b));
    pil_vec2 dd = v_sub(pil_anchor_vel(e->a), pil_anchor_vel(e->b));
    double wr = a_omega(e->a);
    double wb = a_omega(e->b);
    double thr = e->a.body ? e->a.body->angle : 0;
    double thb = e->b.body ? e->b.body->angle : 0;
    double sg = e->sigma, r = e->radius, R = e->radius2, D = R + sg * r;
    double L2 = v_dot(d, d);
    if (L2 < 1e-18) {
      L2 = 1e-18;
    }
    double L = sqrt(L2);
    pil_vec2 n = v_mul(d, 1.0 / L);
    // row 0: centre distance  |C - O| = R + \sigma r
    double ndd = v_dot(n, dd);
    rows[0].phi = L - D;
    rows[0].vel = ndd;
    rows[0].ja[0] = n.x;
    rows[0].ja[1] = n.y;
    rows[0].ja[2] = v_cross(rr, n);
    rows[0].jb[0] = -n.x;
    rows[0].jb[1] = -n.y;
    rows[0].jb[2] = -v_cross(rb, n);
    finish(e, &rows[0],
           v_dot(n, v_sub(v_mul(rr, wr * wr), v_mul(rb, wb * wb))) -
               (v_dot(dd, dd) - ndd * ndd) / L);
    // row 1: no slip.  Phi = D*psi_rel + sigma*r*(th_b - th_r) - c, psi_rel =
    // angle of d in track frame, unwrapped about reference angle (ball never
    // travels more than pi around track)
    pil_vec2 dl = v_rot(d, -thb);
    double raw = atan2(dl.y, dl.x) - e->ref;
    while (raw > 3.14159265358979323846) {
      raw -= 6.28318530717958647692;
    }
    while (raw < -3.14159265358979323846) {
      raw += 6.28318530717958647692;
    }
    double psi_rel = e->ref + raw;
    double psid = v_cross(d, dd) / L2;
    pil_vec2 g = v_mul(v_perp(d), D / L2);
    rows[1].phi = D * psi_rel + sg * r * (thb - thr) - e->phase;
    rows[1].vel = D * psid - sg * r * wr - R * wb;
    rows[1].ja[0] = g.x;
    rows[1].ja[1] = g.y;
    rows[1].ja[2] = v_cross(rr, g) - sg * r;
    rows[1].jb[0] = -g.x;
    rows[1].jb[1] = -g.y;
    rows[1].jb[2] = -v_cross(rb, g) - R;
    finish(e, &rows[1],
           D * (wr * wr * v_cross(d, rr) / L2 - wb * wb * v_cross(d, rb) / L2 +
                2.0 * psid * v_dot(d, dd) / L2));

    return 2;
  }

  case EL_PULLEY: {
    if (e->rope_off && !force_active) {
      return 0;
    }

    pil_vec2 ra = pil_anchor_r(e->a);
    pil_vec2 rb = pil_anchor_r(e->b);
    pil_vec2 da = v_sub(pil_anchor_world(e->a), e->p1);
    pil_vec2 db = v_sub(pil_anchor_world(e->b), e->p2);
    pil_vec2 va = pil_anchor_vel(e->a);
    pil_vec2 vb = pil_anchor_vel(e->b);
    double la = sqrt(v_dot(da, da));
    double lb = sqrt(v_dot(db, db));
    pil_vec2 na = la > 1e-9 ? v_mul(da, 1.0 / la) : pil_v2(0, 1);
    pil_vec2 nb = lb > 1e-9 ? v_mul(db, 1.0 / lb) : pil_v2(0, 1);
    double nva = v_dot(na, va), nvb = v_dot(nb, vb);
    double ca = la > 1e-9 ? (v_dot(va, va) - nva * nva) / la : 0.0;
    double cb = lb > 1e-9 ? (v_dot(vb, vb) - nvb * nvb) / lb : 0.0;
    double wa2 = a_omega(e->a) * a_omega(e->a);
    double wb2 = a_omega(e->b) * a_omega(e->b);

    rows[0].phi = la + lb - e->length;
    rows[0].vel = nva + nvb;
    rows[0].ja[0] = na.x;
    rows[0].ja[1] = na.y;
    rows[0].ja[2] = v_cross(ra, na);
    rows[0].jb[0] = nb.x;
    rows[0].jb[1] = nb.y;
    rows[0].jb[2] = v_cross(rb, nb);
    finish(e, &rows[0], wa2 * v_dot(na, ra) + wb2 * v_dot(nb, rb) - ca - cb);

    return 1;
  }

  case EL_SLOT: {
    pil_vec2 ra = pil_anchor_r(e->a);
    pil_vec2 r0 = pil_anchor_r(e->b);
    double wB = a_omega(e->b);
    double wa = a_omega(e->a);
    pil_vec2 u = v_rot(e->dir, e->b.body ? e->b.body->angle : 0.0);
    pil_vec2 n = v_perp(u);
    pil_vec2 d = v_sub(pil_anchor_world(e->a), pil_anchor_world(e->b));
    pil_vec2 dd = v_sub(pil_anchor_vel(e->a), pil_anchor_vel(e->b));
    double ud = v_dot(u, d);

    rows[0].phi = v_dot(n, d);
    rows[0].vel = -wB * ud + v_dot(n, dd);
    rows[0].ja[0] = n.x;
    rows[0].ja[1] = n.y;
    rows[0].ja[2] = v_cross(ra, n);
    rows[0].jb[0] = -n.x;
    rows[0].jb[1] = -n.y;
    rows[0].jb[2] = -v_cross(r0, n) - ud;
    finish(e, &rows[0],
           wB * wB * v_dot(n, d) + 2 * wB * v_dot(u, dd) +
               wa * wa * v_dot(n, ra) - wB * wB * v_dot(n, r0));

    return 1;
  }

  case EL_FIXED_ANGLE: {
    double ta = e->a.body ? e->a.body->angle : 0;
    double tb = e->b.body ? e->b.body->angle : 0;
    rows[0].phi = ta - tb - e->phase;
    rows[0].vel = a_omega(e->a) - a_omega(e->b);
    rows[0].ja[2] = 1;
    rows[0].jb[2] = -1;
    finish(e, &rows[0], 0.0);

    return 1;
  }

  case EL_LINK: {
    rows[0].phi =
        coord_q(e->a.body, e->ca) - e->k * coord_q(e->b.body, e->cb) - e->phase;
    rows[0].vel = coord_v(e->a.body, e->ca) - e->k * coord_v(e->b.body, e->cb);
    rows[0].ja[e->ca] = 1;
    rows[0].jb[e->cb] = -e->k;
    finish(e, &rows[0], 0.0);

    return 1;
  }

  case EL_ROT_DRIVE: {
    rows[0].phi = e->a.body->angle - e->phase - e->omega * t;
    rows[0].vel = e->a.body->omega - e->omega;
    rows[0].ja[2] = 1;
    finish(e, &rows[0], 0.0);

    return 1;
  }

  case EL_ROLLING: {
    pil_body *b = e->a.body;
    pil_vec2 rC = pil_anchor_r(e->a);
    pil_vec2 C = pil_anchor_world(e->a);
    pil_vec2 vC = pil_anchor_vel(e->a);
    pil_vec2 u = e->dir, n = v_perp(u);
    pil_vec2 dC = v_sub(C, e->b.local);

    double w2 = b->omega * b->omega;
    double R = e->radius;

    rows[0].phi = v_dot(n, dC) - R;
    rows[0].vel = v_dot(n, vC);
    rows[0].ja[0] = n.x;
    rows[0].ja[1] = n.y;
    rows[0].ja[2] = v_cross(rC, n);
    finish(e, &rows[0], w2 * v_dot(n, rC));

    rows[1].phi = v_dot(u, dC) + R * b->angle - e->phase;
    rows[1].vel = v_dot(u, vC) + R * b->omega;
    rows[1].ja[0] = u.x;
    rows[1].ja[1] = u.y;
    rows[1].ja[2] = v_cross(rC, u) + R;
    finish(e, &rows[1], w2 * v_dot(u, rC));

    return 2;
  }

  default:
    return 0;
  }
}
