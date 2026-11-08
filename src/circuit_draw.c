/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* circuit_draw.c - textbook circuit symbols
 *
 * Each element is drawn along straight line between its two nodes: wire, symbol
 * centred on line, wire
 * Everything is in screen pixels so symbols keep their size when figure is
 * zoomed */

#include "pil_circuit_draw.h"
#include <cairo.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static pil_vec2 node_px(const pil_circuit *c, const pil_view *v, pil_node n) {
  return pil_view_to_screen(v, pil_v2(c->nodes[n].x, c->nodes[n].y));
}

static void ink(cairo_t *cr, const pil_style *st) {
  cairo_set_source_rgb(cr, st->ink[0], st->ink[1], st->ink[2]);
}

static void line(cairo_t *cr, double x0, double y0, double x1, double y1) {
  cairo_move_to(cr, x0, y0);
  cairo_line_to(cr, x1, y1);
  cairo_stroke(cr);
}

// Symbol drawn in a local frame: origin at centre, x along a -> b
static void symbol(cairo_t *cr, const pil_celement *e, double len,
                   const pil_style *st, double t) {
  const double h = len / 2;
  cairo_set_line_width(cr, st->line_w);

  switch (e->kind) {
  case CE_RESISTOR: {
    const int n = 6; // zigzag: n peaks between two short leads
    const double amp = 8;
    const double lead = len * 0.1;
    const double span = len - 2 * lead;

    cairo_move_to(cr, -h, 0);
    cairo_line_to(cr, -h + lead, 0);
    for (int k = 0; k < n; k++) {
      cairo_line_to(cr, -h + lead + (k + 0.5) * span / n,
                    (k % 2 ? 1 : -1) * amp);
    }

    cairo_line_to(cr, h - lead, 0);
    cairo_line_to(cr, h, 0);
    cairo_stroke(cr);

    break;
  }
  case CE_SWITCH: {
    int closed = t >= e->t_change ? !e->closed_before : e->closed_before;

    line(cr, -h, 0, -h + 6, 0);
    line(cr, h - 6, 0, h, 0);
    cairo_arc(cr, -h + 6, 0, 2.5, 0, 2 * M_PI);
    cairo_stroke(cr);
    cairo_arc(cr, h - 6, 0, 2.5, 0, 2 * M_PI);
    cairo_stroke(cr);

    double lx = h - 8.5;
    const double ly = 0;

    if (closed) {
      line(cr, -h + 8.5, 0, lx, ly);
    } else {
      line(cr, -h + 8.5, 0, lx - 2, -16);
    }

    break;
  }
  case CE_CAPACITOR: {
    const double gap = 4;
    const double pl = 14;

    line(cr, -h, 0, -gap, 0);
    line(cr, gap, 0, h, 0);
    cairo_set_line_width(cr, st->line_w * 1.8);
    line(cr, -gap, -pl, -gap, pl);
    line(cr, gap, -pl, gap, pl);

    break;
  }
  case CE_MOTOR: {
    const double r = 15;

    line(cr, -h, 0, -r, 0);
    line(cr, r, 0, h, 0);
    cairo_arc(cr, 0, 0, r, 0, 2 * M_PI);
    cairo_stroke(cr);

    break; // letter M is drawn upright by caller
  }
  case CE_VCOIL:
  case CE_INDUCTOR: {
    const int n = 4;
    const double lead = len * 0.12;
    const double d = (len - 2 * lead) / n;

    line(cr, -h, 0, -h + lead, 0);
    line(cr, h - lead, 0, h, 0);
    for (int k = 0; k < n; k++) {
      cairo_new_sub_path(cr);
      cairo_arc(cr, -h + lead + (k + 0.5) * d, 0, d / 2, M_PI, 2 * M_PI);
    }

    cairo_stroke(cr);

    break;
  }
  case CE_VSOURCE: {
    if (e->ac) {
      const double r = 15;

      line(cr, -h, 0, -r, 0);
      line(cr, r, 0, h, 0);
      cairo_arc(cr, 0, 0, r, 0, 2 * M_PI);
      cairo_stroke(cr);
      cairo_move_to(cr, -9, 0);
      cairo_curve_to(cr, -5, -12, -2, -12, 0, 0);
      cairo_curve_to(cr, 2, 12, 5, 12, 9, 0);
      cairo_stroke(cr);
    } else {
      const double gap = 4;

      line(cr, -h, 0, -gap, 0);
      line(cr, gap, 0, h, 0);
      cairo_set_line_width(cr, st->line_w * 1.5);
      line(cr, -gap, -15, -gap, 15);
      line(cr, gap, -8, gap, 8);
    }

    break;
  }
  case CE_ISOURCE: {
    const double r = 15;

    line(cr, -h, 0, -r, 0);
    line(cr, r, 0, h, 0);
    cairo_arc(cr, 0, 0, r, 0, 2 * M_PI);
    cairo_stroke(cr);
    line(cr, -9, 0, 7, 0);
    cairo_move_to(cr, 8, 0);
    cairo_line_to(cr, 1, -4.5);
    cairo_line_to(cr, 1, 4.5);
    cairo_close_path(cr);
    cairo_fill(cr);

    break;
  }
  case CE_DIODE: {
    const double w = 9;
    const double tip = 8;

    line(cr, -h, 0, -tip, 0);
    line(cr, tip, 0, h, 0);
    cairo_move_to(cr, -tip, -w);
    cairo_line_to(cr, -tip, w);
    cairo_line_to(cr, tip, 0);
    cairo_close_path(cr);
    cairo_stroke(cr);
    cairo_set_line_width(cr, st->line_w * 1.5);
    line(cr, tip, -w, tip, w);

    break;
  }
  }
}

static double symbol_length(const pil_celement *e) {
  switch (e->kind) {
  case CE_RESISTOR:
    return 58;
  case CE_INDUCTOR:
  case CE_VCOIL:
    return 64;
  case CE_MOTOR:
    return 34;
  case CE_CAPACITOR:
    return 24;
  case CE_VSOURCE:
    return e->ac ? 34 : 24;
  case CE_ISOURCE:
    return 34;
  case CE_DIODE:
    return 34;
  case CE_SWITCH:
    return 52;
  }

  return 40;
}

// One flow mark: an arrow head with its tip at (x, y) pointing along (dx, dy),
// followed by three dots that shrink away behind it
// Light fill, ink rim, so mark reads over both wire and paper
static void flow_mark(cairo_t *cr, double x, double y, double dx, double dy,
                      const pil_style *st) {
  const double r = st->current_dot_r;
  cairo_save(cr);
  cairo_set_line_width(cr, 1.0);
  cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);

  double nx = -dy, ny = dx;
  double hl = 2.8 * r, hw = 1.5 * r; // head length, half width
  cairo_move_to(cr, x, y);
  cairo_line_to(cr, x - dx * hl + nx * hw, y - dy * hl + ny * hw);
  cairo_line_to(cr, x - dx * hl * 0.75, y - dy * hl * 0.75);
  cairo_line_to(cr, x - dx * hl - nx * hw, y - dy * hl - ny * hw);
  cairo_close_path(cr);
  cairo_set_source_rgb(cr, st->current_color[0], st->current_color[1],
                       st->current_color[2]);
  cairo_fill_preserve(cr);
  cairo_set_source_rgb(cr, st->ink[0], st->ink[1], st->ink[2]);
  cairo_stroke(cr);

  for (int j = 0; j < 3; j++) {
    double back = hl + (1.2 + 2.3 * j) * r;
    double rj = r * (1.0 - 0.18 * j);
    cairo_arc(cr, x - dx * back, y - dy * back, rj, 0, 2 * M_PI);
    cairo_set_source_rgb(cr, st->current_color[0], st->current_color[1],
                         st->current_color[2]);
    cairo_fill_preserve(cr);
    cairo_set_source_rgb(cr, st->ink[0], st->ink[1], st->ink[2]);
    cairo_stroke(cr);
  }

  cairo_restore(cr);
}

// Current along branch a -> b. phase_px advances with current (it is
// sign-carrying integral of i / peak), so marks stream towards b for i > 0 and
// towards a for i < 0, faster for bigger currents. skip_mid keeps symbol in
// middle clear
static void draw_current_flow(cairo_t *cr, pil_vec2 a, pil_vec2 b,
                              double phase_px, double dir, double skip_mid,
                              const pil_style *st) {
  double dx = b.x - a.x;
  double dy = b.y - a.y;
  double L = hypot(dx, dy);
  if (L < 14) {
    return;
  }

  double ux = dx / L;
  double uy = dy / L;
  double mid_lo = (L - skip_mid) * 0.5;
  double mid_hi = mid_lo + skip_mid;
  const double period = 13.0 * st->current_dot_r;
  double off = fmod(phase_px, period);
  if (off < 0) {
    off += period;
  }

  // Mark occupies [tip - extent, tip] along flow
  // stay inside branch and outside symbol
  double extent = 10.6 * st->current_dot_r;
  for (double s = off; s < L + period; s += period) {
    double lo = dir > 0 ? s - extent : s;
    double hi = dir > 0 ? s : s + extent;
    if (lo < 1 || hi > L - 1) {
      continue;
    }

    if (skip_mid > 0 && hi > mid_lo && lo < mid_hi) {
      continue;
    }

    flow_mark(cr, a.x + ux * s, a.y + uy * s, ux * dir, uy * dir, st);
  }
}

static void voltage_text(char *buf, size_t n, double v) {
  if (fabs(v) < 5e-4) {
    snprintf(buf, n, "0 V");
  } else if (fabs(v) < 1.0) {
    snprintf(buf, n, "%.0f mV", v * 1e3);
  } else {
    snprintf(buf, n, "%.3g V", v);
  }
}

static void draw_ground(cairo_t *cr, pil_vec2 p, const pil_style *st) {
  cairo_set_line_width(cr, st->line_w);
  line(cr, p.x, p.y, p.x, p.y + 10);
  line(cr, p.x - 12, p.y + 10, p.x + 12, p.y + 10);
  line(cr, p.x - 8, p.y + 15, p.x + 8, p.y + 15);
  line(cr, p.x - 4, p.y + 20, p.x + 4, p.y + 20);
}

void pil_circuit_draw(cairo_t *cr, const pil_circuit *c, const pil_view *v,
                      const pil_style *st) {
  cairo_save(cr);
  ink(cr, st);
  cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
  cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
  cairo_set_line_width(cr, st->line_w);

  int *deg = calloc((size_t)c->nn + 1, sizeof *deg);
  if (!deg) {
    abort();
  }

  // marks only for currents that matter next to the strongest one
  double i_min = 2e-3 * c->peak_i;
  int flowing =
      st->show_current && c->peak_i > 1e-12 && st->current_dot_speed > 0;

  for (int k = 0; k < c->nw; k++) {
    pil_vec2 a = node_px(c, v, c->wires[k].a);
    pil_vec2 b = node_px(c, v, c->wires[k].b);
    line(cr, a.x, a.y, b.x, b.y);
    deg[c->wires[k].a]++;
    deg[c->wires[k].b]++;

    if (flowing && fabs(c->wires[k].i) > i_min) {
      draw_current_flow(cr, a, b,
                        c->wires[k].anim_phase * st->current_dot_speed,
                        c->wires[k].i > 0 ? 1 : -1, 0, st);
    }
  }

  for (int k = 0; k < c->ne; k++) {
    const pil_celement *e = c->els[k];
    pil_vec2 a = node_px(c, v, e->a);
    pil_vec2 b = node_px(c, v, e->b);
    deg[e->a]++;
    deg[e->b]++;
    double L = hypot(b.x - a.x, b.y - a.y);
    if (L < 1e-6) {
      continue;
    }

    double ux = (b.x - a.x) / L;
    double uy = (b.y - a.y) / L;
    double sl = fmin(symbol_length(e), 0.8 * L);
    pil_vec2 m = {0.5 * (a.x + b.x), 0.5 * (a.y + b.y)};
    line(cr, a.x, a.y, m.x - ux * sl / 2, m.y - uy * sl / 2);
    line(cr, m.x + ux * sl / 2, m.y + uy * sl / 2, b.x, b.y);

    if (!e->hidden) {
      cairo_save(cr);
      cairo_translate(cr, m.x, m.y);
      cairo_rotate(cr, atan2(uy, ux));
      symbol(cr, e, sl, st, c->world->t);
      cairo_restore(cr);
    }

    if (flowing && fabs(e->i) > i_min) {
      draw_current_flow(cr, a, b, e->anim_phase * st->current_dot_speed,
                        e->i > 0 ? 1 : -1, sl, st);
    }

    if (e->kind == CE_MOTOR && !e->hidden) {
      pil_draw_text(cr, m, "M", 0.5, st);
    }

    if (e->label[0] && !e->hidden) {
      double dx = e->label_dx;
      double dy = e->label_dy;
      if (dx == 0 && dy == 18 && fabs(uy) > fabs(ux)) {
        dx = 26; // default for a vertical element: beside it
        dy = 0;
      }

      pil_draw_text(cr, pil_v2(m.x + dx, m.y - dy), e->label,
                    dx > 8 ? 0.0 : (dx < -8 ? 1.0 : 0.5), st);
    }
  }

  for (int n = 0; n < c->nn; n++) {
    if (deg[n] >= 3) {
      pil_vec2 p = node_px(c, v, n);
      cairo_arc(cr, p.x, p.y, 3.2, 0, 2 * M_PI);
      cairo_fill(cr);
    }
  }

  // node voltages at joints, once circuit is running. One label per electrical
  // node (wires merge nodes) where at least two elements meet
  if (st->show_voltages && c->peak_i > 1e-12 && c->vn) {
    int *terms = calloc((size_t)c->nn, sizeof *terms);
    if (terms) {
      for (int k = 0; k < c->ne; k++) {
        terms[pil_cnode_find((pil_circuit *)c, c->els[k]->a)]++;
        terms[pil_cnode_find((pil_circuit *)c, c->els[k]->b)]++;
      }

      int gr =
          c->ground >= 0 ? pil_cnode_find((pil_circuit *)c, c->ground) : -1;
      for (int n = 0; n < c->nn; n++) {
        if (pil_cnode_find((pil_circuit *)c, n) != n || n == gr ||
            terms[n] < 2) {
          continue;
        }

        char buf[24];
        voltage_text(buf, sizeof buf, c->vn[n]);
        pil_vec2 p = node_px(c, v, n);
        pil_draw_text(cr, pil_v2(p.x + 6, p.y - 17), buf, 0.0, st);
      }

      free(terms);
    }
  }

  if (c->ground >= 0) {
    draw_ground(cr, node_px(c, v, c->ground), st);
  }

  free(deg);
  cairo_restore(cr);
}
