/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* circuit_draw.c - textbook circuit symbols
 *
 * Each element is drawn along straight line between its two nodes: wire, symbol
 * centred on line, wire
 * Everything is in screen pixels so symbols keep their size when figure is
 * zoomed
 */

#include "pil_circuit_draw.h"
#include <cairo.h>
#include <math.h>
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
    return 64;
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

  for (int k = 0; k < c->nw; k++) {
    pil_vec2 a = node_px(c, v, c->wires[k].a), b = node_px(c, v, c->wires[k].b);
    line(cr, a.x, a.y, b.x, b.y);
    deg[c->wires[k].a]++;
    deg[c->wires[k].b]++;
  }

  for (int k = 0; k < c->ne; k++) {
    const pil_celement *e = c->els[k];
    pil_vec2 a = node_px(c, v, e->a), b = node_px(c, v, e->b);
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

  if (c->ground >= 0) {
    draw_ground(cr, node_px(c, v, c->ground), st);
  }

  free(deg);
  cairo_restore(cr);
}
