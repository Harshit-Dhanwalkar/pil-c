/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* draw.c - textbook-style Cairo renderer */

#include "pil/pil_draw.h"
#include "pil_internal.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define TWO_PI 6.28318530717958647692

// style & view
pil_style pil_style_textbook(void) {
  pil_style s;
  memset(&s, 0, sizeof s);

  s.line_w = 1.5;
  s.thin_w = 1.1;
  s.body_fill = 0.82;
  s.font_px = 18;
  s.font = "Latin Modern Roman";
  s.math_font = "Latin Modern Math";
  s.ink[0] = s.ink[1] = s.ink[2] = 0.05;
  s.paper[0] = s.paper[1] = s.paper[2] = 1.0;
  s.draw_background = 1;
  s.show_com = 1;
  s.spring_coils = 7;
  s.spring_w_px = 10.0;
  s.hinge_r_px = 3.6;

  return s;
}

void pil_view_fit(pil_view *v, double x0, double y0, double x1, double y1,
                  double W, double H, double margin) {
  double sx = (W - 2 * margin) / (x1 - x0);
  double sy = (H - 2 * margin) / (y1 - y0);
  v->scale = sx < sy ? sx : sy;
  v->cx = 0.5 * (x0 + x1);
  v->cy = 0.5 * (y0 + y1);
  v->width = W;
  v->height = H;
}

pil_vec2 pil_view_to_screen(const pil_view *v, pil_vec2 p) {
  return pil_v2(0.5 * v->width + (p.x - v->cx) * v->scale,
                0.5 * v->height - (p.y - v->cy) * v->scale);
}

pil_vec2 pil_view_to_world(const pil_view *v, pil_vec2 s) {
  return pil_v2(v->cx + (s.x - 0.5 * v->width) / v->scale,
                v->cy - (s.y - 0.5 * v->height) / v->scale);
}

static void ink(cairo_t *cr, const pil_style *st) {
  cairo_set_source_rgb(cr, st->ink[0], st->ink[1], st->ink[2]);
}
static void paper(cairo_t *cr, const pil_style *st) {
  cairo_set_source_rgb(cr, st->paper[0], st->paper[1], st->paper[2]);
}

static pil_vec2 S(const pil_view *v, pil_vec2 w) {
  return pil_view_to_screen(v, w);
}

// text with sub/superscripts
static int utf8_len(unsigned char c) {
  return c < 0x80         ? 1
         : (c >> 5) == 6  ? 2
         : (c >> 4) == 14 ? 3
         : (c >> 3) == 30 ? 4
                          : 1;
}

static unsigned utf8_decode(const char *s, int n) {
  unsigned char c = (unsigned char)s[0];
  if (n == 1) {
    return c;
  }

  unsigned cp = n == 2 ? (c & 0x1f) : n == 3 ? (c & 0x0f) : (c & 0x07);
  for (int i = 1; i < n; i++) {
    cp = (cp << 6) | ((unsigned char)s[i] & 0x3f);
  }

  return cp;
}

static int utf8_encode(unsigned cp, char *o) {
  if (cp < 0x80) {
    o[0] = (char)cp;
    return 1;
  }
  if (cp < 0x800) {
    o[0] = (char)(0xc0 | (cp >> 6));
    o[1] = (char)(0x80 | (cp & 0x3f));
    return 2;
  }
  if (cp < 0x10000) {
    o[0] = (char)(0xe0 | (cp >> 12));
    o[1] = (char)(0x80 | ((cp >> 6) & 0x3f));
    o[2] = (char)(0x80 | (cp & 0x3f));
    return 3;
  }
  o[0] = (char)(0xf0 | (cp >> 18));
  o[1] = (char)(0x80 | ((cp >> 12) & 0x3f));
  o[2] = (char)(0x80 | ((cp >> 6) & 0x3f));
  o[3] = (char)(0x80 | (cp & 0x3f));

  return 4;
}

// Greek letters and script l are not in Latin Modern Roman: draw them from
// maths font, mapped to Unicode "mathematical italic" block
static int is_math_glyph(unsigned cp) {
  return cp == 0x2113 || (cp >= 0x391 && cp <= 0x3a9) ||
         (cp >= 0x3b1 && cp <= 0x3c9);
}

static unsigned math_italic(unsigned cp) {
  if (cp >= 0x391 && cp <= 0x3a9) {
    return 0x1D6E2 + (cp - 0x391);
  }
  if (cp >= 0x3b1 && cp <= 0x3c9) {
    return 0x1D6FC + (cp - 0x3b1);
  }

  return cp;
}

static double text_piece(cairo_t *cr, const char *family,
                         cairo_font_slant_t slant, const char *s, int n,
                         double x, double y, double size, int draw) {
  char buf[128];
  if (n <= 0) {
    return 0;
  }
  if (n > 127) {
    n = 127;
  }

  memcpy(buf, s, (size_t)n);
  buf[n] = 0;
  cairo_select_font_face(cr, family, slant, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, size);
  cairo_text_extents_t ext;
  cairo_text_extents(cr, buf, &ext);
  if (draw) {
    cairo_move_to(cr, x, y);
    cairo_show_text(cr, buf);
  }

  return ext.x_advance;
}

// one chunk (no _ or ^): splits into Latin runs and maths glyphs
static double text_chunk(cairo_t *cr, const pil_style *st, const char *s, int n,
                         double x, double y, double size, int draw) {
  double x0 = x;
  int i = 0;
  while (i < n) {
    int j = i;
    while (j < n) {
      int l = utf8_len((unsigned char)s[j]);
      if (is_math_glyph(utf8_decode(s + j, l))) {
        break;
      }
      j += l;
    }

    if (j > i)
      x += text_piece(cr, st->font, CAIRO_FONT_SLANT_ITALIC, s + i, j - i, x, y,
                      size, draw);
    if (j < n) {
      int l = utf8_len((unsigned char)s[j]);
      char g[8];
      int gl = utf8_encode(math_italic(utf8_decode(s + j, l)), g);
      x += text_piece(cr, st->math_font, CAIRO_FONT_SLANT_NORMAL, g, gl, x, y,
                      size, draw);
      j += l;
    }
    i = j;
  }

  return x - x0;
}

static double text_run(cairo_t *cr, const pil_style *st, const char *s,
                       double x, double y, int draw) {
  double x0 = x, fs = st->font_px;
  while (*s) {
    if (*s == '_' || *s == '^') {
      int sup = (*s == '^');
      s++;
      const char *a = s;
      int n;
      if (*s == '{') {
        a = ++s;
        while (*s && *s != '}') {
          s++;
        }
        n = (int)(s - a);
        if (*s == '}') {
          s++;
        }
      } else {
        n = utf8_len((unsigned char)*s);
        s += n;
      }
      x += text_chunk(cr, st, a, n, x, y + (sup ? -0.38 : 0.24) * fs, 0.7 * fs,
                      draw);
    } else {
      const char *a = s;
      while (*s && *s != '_' && *s != '^') {
        s++;
      }
      x += text_chunk(cr, st, a, (int)(s - a), x, y, fs, draw);
    }
  }

  return x - x0;
}

void pil_draw_text(cairo_t *cr, pil_vec2 p, const char *text, double align,
                   const pil_style *st) {
  if (!text || !*text) {
    return;
  }

  cairo_save(cr);
  double w = text_run(cr, st, text, 0, 0, 0);
  ink(cr, st);
  text_run(cr, st, text, p.x - align * w, p.y + 0.33 * st->font_px, 1);
  cairo_restore(cr);
}

static double text_width(cairo_t *cr, const pil_style *st, const char *t) {
  cairo_save(cr);
  double w = text_run(cr, st, t, 0, 0, 0);
  cairo_restore(cr);

  return w;
}

/* symbols*/
void pil_draw_hinge(cairo_t *cr, pil_vec2 p, const pil_style *st) {
  cairo_save(cr);
  cairo_new_path(cr);
  cairo_arc(cr, p.x, p.y, st->hinge_r_px, 0, TWO_PI);
  paper(cr, st);
  cairo_fill_preserve(cr);
  ink(cr, st);
  cairo_set_line_width(cr, st->thin_w * 1.2);
  cairo_stroke(cr);
  cairo_restore(cr);
}

void pil_draw_com(cairo_t *cr, pil_vec2 p, const pil_style *st) {
  const double r = 5.5;
  cairo_save(cr);
  cairo_new_path(cr);
  cairo_arc(cr, p.x, p.y, r, 0, TWO_PI);
  paper(cr, st);
  cairo_fill(cr);
  ink(cr, st);
  cairo_move_to(cr, p.x, p.y);
  cairo_arc(cr, p.x, p.y, r, 0, 0.25 * TWO_PI);
  cairo_close_path(cr);
  cairo_fill(cr);
  cairo_move_to(cr, p.x, p.y);
  cairo_arc(cr, p.x, p.y, r, 0.5 * TWO_PI, 0.75 * TWO_PI);
  cairo_close_path(cr);
  cairo_fill(cr);
  cairo_new_path(cr);
  cairo_arc(cr, p.x, p.y, r, 0, TWO_PI);
  cairo_set_line_width(cr, st->thin_w * 1.1);
  cairo_stroke(cr);
  cairo_restore(cr);
}

void pil_draw_arrow(cairo_t *cr, pil_vec2 a, pil_vec2 b, double head,
                    const pil_style *st) {
  double dx = b.x - a.x;
  double dy = b.y - a.y;
  double L = sqrt(dx * dx + dy * dy);
  if (L < 1e-6) {
    return;
  }

  dx /= L;
  dy /= L;
  if (head > 0.6 * L) {
    head = 0.6 * L;
  }

  cairo_save(cr);
  ink(cr, st);
  cairo_set_line_width(cr, st->thin_w * 1.15);
  cairo_set_line_cap(cr, CAIRO_LINE_CAP_BUTT);
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, b.x - dx * head * 0.8, b.y - dy * head * 0.8);
  cairo_stroke(cr);
  cairo_move_to(cr, b.x, b.y);
  cairo_line_to(cr, b.x - dx * head - dy * head * 0.3,
                b.y - dy * head + dx * head * 0.3);
  cairo_line_to(cr, b.x - dx * head * 0.75, b.y - dy * head * 0.75);
  cairo_line_to(cr, b.x - dx * head + dy * head * 0.3,
                b.y - dy * head - dx * head * 0.3);
  cairo_close_path(cr);
  cairo_fill(cr);
  cairo_restore(cr);
}

void pil_draw_spring(cairo_t *cr, pil_vec2 a, pil_vec2 b, int coils, double hw,
                     const pil_style *st) {
  double dx = b.x - a.x;
  double dy = b.y - a.y;
  double L = sqrt(dx * dx + dy * dy);
  if (L < 2) {
    return;
  }

  dx /= L;
  dy /= L;
  double px = -dy, py = dx;
  double lead = 0.17 * L;
  double Lc = L - 2 * lead;
  double bb = 0.85 * hw;
  int per = 28, N = coils * per;
  cairo_save(cr);
  ink(cr, st);
  cairo_set_line_width(cr, st->thin_w * 1.25);
  cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
  cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, a.x + dx * lead, a.y + dy * lead);
  for (int i = 1; i <= N; i++) {
    double t = (double)i / N, th = TWO_PI * coils * t;
    double along = lead + Lc * t + bb * (1.0 - cos(th));
    double lat = hw * sin(th);
    cairo_line_to(cr, a.x + dx * along + px * lat, a.y + dy * along + py * lat);
  }
  cairo_line_to(cr, b.x, b.y);
  cairo_stroke(cr);
  cairo_restore(cr);
}

void pil_draw_dashpot(cairo_t *cr, pil_vec2 a, pil_vec2 b, double nominal,
                      double hw, const pil_style *st) {
  double dx = b.x - a.x;
  double dy = b.y - a.y;
  double L = sqrt(dx * dx + dy * dy);
  if (L < 2) {
    return;
  }

  dx /= L;
  dy /= L;
  double px = -dy, py = dx;
  if (nominal < 20) {
    nominal = L;
  }

  double s0 = 0.18 * nominal, s1 = s0 + 0.5 * nominal;
  double sp = L - 0.32 * nominal;
  if (sp < s0 + 3) {
    sp = s0 + 3;
  }

  double w = hw * 0.95;
#define P(s, l) a.x + dx *(s) + px *(l), a.y + dy *(s) + py *(l)
  cairo_save(cr);
  ink(cr, st);
  cairo_set_line_width(cr, st->thin_w * 1.25);
  cairo_set_line_join(cr, CAIRO_LINE_JOIN_MITER);
  cairo_set_line_cap(cr, CAIRO_LINE_CAP_BUTT);
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, P(s0, 0));
  cairo_stroke(cr); /* lead to cylinder */
  cairo_move_to(cr, P(s1, -w));
  cairo_line_to(cr, P(s0, -w));
  cairo_line_to(cr, P(s0, w));
  cairo_line_to(cr, P(s1, w));
  cairo_stroke(cr); /* U-shaped cylinder */
  cairo_move_to(cr, P(sp, -w * 0.8));
  cairo_line_to(cr, P(sp, w * 0.8));
  cairo_stroke(cr); /* piston */
  cairo_move_to(cr, P(sp, 0));
  cairo_line_to(cr, b.x, b.y);
  cairo_stroke(cr); /* rod */
  cairo_restore(cr);
#undef P
}

void pil_draw_hatch(cairo_t *cr, pil_vec2 a, pil_vec2 b, const pil_style *st) {
  double dx = b.x - a.x;
  double dy = b.y - a.y;
  double L = sqrt(dx * dx + dy * dy);
  if (L < 1) {
    return;
  }

  dx /= L;
  dy /= L;
  double rx = -dy, ry = dx; // right of travel on screen (y down)
  cairo_save(cr);
  ink(cr, st);
  cairo_set_line_width(cr, st->line_w);
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, b.x, b.y);
  cairo_stroke(cr);
  cairo_set_line_width(cr, st->thin_w);
  for (double s = 0; s <= L; s += 9.0) {
    double x = a.x + dx * s;
    double y = a.y + dy * s;
    cairo_move_to(cr, x, y);
    cairo_line_to(cr, x + (rx - dx) * 7, y + (ry - dy) * 7);
    cairo_stroke(cr);
  }
  cairo_restore(cr);
}

/* bodies */
static void body_path(cairo_t *cr, const pil_body *b, const pil_view *v) {
  const pil_shape *s = &b->shape;
  pil_vec2 c0 = v_add(b->pos, v_rot(s->offset, b->angle));
  cairo_new_path(cr);
  switch (s->kind) {
  case pil_SHAPE_RECT:
    for (int i = 0; i < 4; i++) {
      static const int cx[4] = {-1, 1, 1, -1}, cy[4] = {-1, -1, 1, 1};
      pil_vec2 p =
          S(v, v_add(c0, v_rot(pil_v2(0.5 * cx[i] * s->w, 0.5 * cy[i] * s->h),
                               b->angle)));
      if (i == 0) {
        cairo_move_to(cr, p.x, p.y);
      } else {
        cairo_line_to(cr, p.x, p.y);
      }
    }
    cairo_close_path(cr);
    break;
  case pil_SHAPE_POLY:
    for (int i = 0; i < s->npts; i++) {
      pil_vec2 p = S(v, v_add(c0, v_rot(s->pts[i], b->angle)));
      if (i == 0) {
        cairo_move_to(cr, p.x, p.y);
      } else {
        cairo_line_to(cr, p.x, p.y);
      }
    }
    cairo_close_path(cr);
    break;
  case pil_SHAPE_DISC: {
    pil_vec2 c = S(v, c0);
    cairo_arc(cr, c.x, c.y, s->r * v->scale, 0, TWO_PI);
    break;
  }
  case pil_SHAPE_POINT: {
    pil_vec2 c = S(v, c0);
    double r = s->r * v->scale;
    if (r < 3.5) {
      r = 3.5;
    }
    cairo_arc(cr, c.x, c.y, r, 0, TWO_PI);
    break;
  }
  case pil_SHAPE_ARC_BAND: {
    int N = 72, first = 1;
    for (int i = 0; i <= N; i++) {
      double a = s->a0 + (s->a1 - s->a0) * i / N + b->angle;
      pil_vec2 p =
          S(v, v_add(c0, pil_v2(s->r_out * cos(a), s->r_out * sin(a))));
      if (first) {
        cairo_move_to(cr, p.x, p.y);
        first = 0;
      } else {
        cairo_line_to(cr, p.x, p.y);
      }
    }
    for (int i = N; i >= 0; i--) {
      double a = s->a0 + (s->a1 - s->a0) * i / N + b->angle;
      pil_vec2 p = S(v, v_add(c0, pil_v2(s->r_in * cos(a), s->r_in * sin(a))));
      cairo_line_to(cr, p.x, p.y);
    }
    cairo_close_path(cr);
    break;
  }
  }
}

static void draw_body(cairo_t *cr, const pil_body *b, const pil_view *v,
                      const pil_style *st) {
  if (b->flags & pil_HIDDEN) {
    return;
  }
  cairo_save(cr);
  body_path(cr, b, v);
  if (b->shape.kind == pil_SHAPE_POINT) {
    ink(cr, st);
    cairo_fill(cr);
  } else {
    if (!(b->flags & pil_NO_FILL)) {
      double g = b->fill > 0 ? b->fill : st->body_fill;
      cairo_set_source_rgb(cr, g, g, g);
      cairo_fill_preserve(cr);
    }
    ink(cr, st);
    cairo_set_line_width(cr, st->line_w);
    cairo_set_line_join(cr, CAIRO_LINE_JOIN_MITER);
    cairo_stroke(cr);
  }

  if (b->shape.kind == pil_SHAPE_DISC && (b->flags & pil_SPOKE)) {
    pil_vec2 c0 = v_add(b->pos, v_rot(b->shape.offset, b->angle));
    pil_vec2 c = S(v, c0),
             e = S(v, v_add(c0, v_rot(pil_v2(b->shape.r * 0.9, 0), b->angle)));
    ink(cr, st);
    cairo_set_line_width(cr, st->thin_w);
    cairo_move_to(cr, c.x, c.y);
    cairo_line_to(cr, e.x, e.y);
    cairo_stroke(cr);
  }
  cairo_restore(cr);
}

static double align_for(double dx) {
  return dx < -1 ? 1.0 : (dx > 1 ? 0.0 : 0.5);
}

/* elements */
static void element_label(cairo_t *cr, const pil_element *e, pil_vec2 mid,
                          const pil_style *st) {
  if (!e->label[0]) {
    return;
  }
  pil_vec2 p = pil_v2(mid.x + e->label_dx, mid.y - e->label_dy);
  pil_draw_text(cr, p, e->label, align_for(e->label_dx), st);
}

static void draw_element(cairo_t *cr, const pil_element *e, const pil_view *v,
                         const pil_style *st) {
  if (e->hidden) {
    return;
  }
  switch (e->kind) {
  case EL_SPRING:
  case EL_DAMPER: {
    pil_vec2 a = S(v, pil_anchor_world(e->a)), b = S(v, pil_anchor_world(e->b));
    double hw = e->width_px > 0 ? e->width_px : st->spring_w_px;
    if (e->kind == EL_SPRING) {
      pil_draw_spring(cr, a, b, e->coils > 0 ? e->coils : st->spring_coils, hw,
                      st);
    } else {
      pil_draw_dashpot(cr, a, b, e->length * v->scale, hw, st);
    }
    element_label(cr, e, pil_v2(0.5 * (a.x + b.x), 0.5 * (a.y + b.y)), st);
    break;
  }
  case EL_DISTANCE: {
    pil_vec2 a = S(v, pil_anchor_world(e->a)), b = S(v, pil_anchor_world(e->b));
    cairo_save(cr);
    if (e->flags & pil_ROPE) {
      ink(cr, st);
      cairo_set_line_width(cr, st->thin_w * 1.3);
      cairo_move_to(cr, a.x, a.y);
      cairo_line_to(cr, b.x, b.y);
      cairo_stroke(cr);
    } else {
      cairo_set_line_cap(cr, CAIRO_LINE_CAP_BUTT);
      ink(cr, st);
      cairo_set_line_width(cr, 4.4);
      cairo_move_to(cr, a.x, a.y);
      cairo_line_to(cr, b.x, b.y);
      cairo_stroke(cr);
      paper(cr, st);
      cairo_set_line_width(cr, 2.0);
      cairo_move_to(cr, a.x, a.y);
      cairo_line_to(cr, b.x, b.y);
      cairo_stroke(cr);
    }
    cairo_restore(cr);
    element_label(cr, e, pil_v2(0.5 * (a.x + b.x), 0.5 * (a.y + b.y)), st);
    break;
  }
  case EL_PULLEY: {
    pil_vec2 a = S(v, pil_anchor_world(e->a)), b = S(v, pil_anchor_world(e->b));
    pil_vec2 c1 = S(v, e->p1), c2 = S(v, e->p2);
    cairo_save(cr);
    ink(cr, st);
    cairo_set_line_width(cr, st->thin_w * 1.3);
    cairo_move_to(cr, a.x, a.y);
    cairo_line_to(cr, c1.x, c1.y);
    cairo_line_to(cr, c2.x, c2.y);
    cairo_line_to(cr, b.x, b.y);
    cairo_stroke(cr);
    cairo_restore(cr);
    element_label(cr, e, pil_v2(0.5 * (a.x + c1.x), 0.5 * (a.y + c1.y)), st);
    break;
  }
  default:
    break;
  }
}

static void draw_joint_marks(cairo_t *cr, const pil_element *e,
                             const pil_view *v, const pil_style *st) {
  if (e->hidden) {
    return;
  }
  if (e->kind == EL_REVOLUTE) {
    pil_draw_hinge(cr, S(v, pil_anchor_world(e->a)), st);
  } else if (e->kind == EL_DISTANCE && !(e->flags & pil_ROPE)) {
    pil_draw_hinge(cr, S(v, pil_anchor_world(e->a)), st);
    pil_draw_hinge(cr, S(v, pil_anchor_world(e->b)), st);
  } else if (e->kind == EL_SPRING || e->kind == EL_DAMPER) {
    // small pin at attachment points, as in the figures
    if (e->a.body) {
      pil_draw_hinge(cr, S(v, pil_anchor_world(e->a)), st);
    }
    if (e->b.body) {
      pil_draw_hinge(cr, S(v, pil_anchor_world(e->b)), st);
    }
  }
}

static void draw_terrain(cairo_t *cr, const pil_terrain *t, double time,
                         const pil_view *v, const pil_style *st) {
  cairo_save(cr);
  ink(cr, st);
  cairo_set_line_width(cr, st->line_w);
  cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
  for (double sx = -10; sx <= v->width + 10; sx += 2.0) {
    pil_vec2 w = pil_view_to_world(v, pil_v2(sx, 0));
    pil_vec2 p = S(v, pil_v2(w.x, pil_terrain_height(t, w.x, time)));
    if (sx < 0) {
      cairo_move_to(cr, p.x, p.y);
    } else {
      cairo_line_to(cr, p.x, p.y);
    }
  }

  cairo_stroke(cr);
  cairo_set_line_width(cr, st->thin_w);
  for (double sx = 0; sx <= v->width; sx += 9.0) {
    pil_vec2 w = pil_view_to_world(v, pil_v2(sx, 0));
    pil_vec2 p = S(v, pil_v2(w.x, pil_terrain_height(t, w.x, time)));
    cairo_move_to(cr, p.x, p.y);
    cairo_line_to(cr, p.x - 7, p.y + 8);
    cairo_stroke(cr);
  }
  cairo_restore(cr);
}

/* annotations                                                         */
static void draw_deco(cairo_t *cr, const pil_deco *d, const pil_view *v,
                      const pil_style *st) {
  if (d->hidden) {
    return;
  }
  pil_vec2 A = S(v, pil_anchor_world(d->a)), B = S(v, pil_anchor_world(d->b));
  switch (d->kind) {
  case D_LABEL:
    pil_draw_text(cr, pil_v2(A.x + d->dx, A.y - d->dy), d->text,
                  align_for(d->dx), st);
    break;
  case D_GROUND:
    pil_draw_hatch(cr, A, B, st);
    break;
  case D_GRAVITY:
    pil_draw_arrow(cr, A, pil_v2(A.x, A.y + 38), 10, st);
    pil_draw_text(cr, pil_v2(A.x + 9, A.y + 14), d->text, 0.0, st);
    break;
  case D_LINE:
    cairo_save(cr);
    ink(cr, st);
    cairo_set_line_width(cr, st->thin_w * 0.9);
    if (d->dashed) {
      double dash[4] = {14, 3.5, 2, 3.5};
      cairo_set_dash(cr, dash, 4, 0);
    }
    cairo_move_to(cr, A.x, A.y);
    cairo_line_to(cr, B.x, B.y);
    cairo_stroke(cr);
    cairo_restore(cr);
    break;
  case D_ARROW: {
    pil_vec2 T = pil_v2(A.x + d->dx, A.y - d->dy);
    pil_draw_arrow(cr, A, T, 10, st);
    if (d->text[0]) {
      double L = sqrt(d->dx * d->dx + d->dy * d->dy);
      double ux = L > 0 ? d->dx / L : 1, uy = L > 0 ? -d->dy / L : 0;
      pil_draw_text(cr, pil_v2(T.x + ux * 10, T.y + uy * 10), d->text,
                    ux < -0.3 ? 1.0 : (ux > 0.3 ? 0.0 : 0.5), st);
    }
    break;
  }
  case D_DIM: {
    double dx = B.x - A.x, dy = B.y - A.y, L = sqrt(dx * dx + dy * dy);
    if (L < 4) {
      break;
    }
    dx /= L;
    dy /= L;
    double nx = dy, ny = -dx, o = d->offset, ext = o + (o >= 0 ? 5 : -5);
    pil_vec2 A2 = pil_v2(A.x + nx * o, A.y + ny * o),
             B2 = pil_v2(B.x + nx * o, B.y + ny * o);
    cairo_save(cr);
    ink(cr, st);
    cairo_set_line_width(cr, st->thin_w * 0.85);
    cairo_move_to(cr, A.x + nx * 3 * (o >= 0 ? 1 : -1),
                  A.y + ny * 3 * (o >= 0 ? 1 : -1));
    cairo_line_to(cr, A.x + nx * ext, A.y + ny * ext);
    cairo_stroke(cr);
    cairo_move_to(cr, B.x + nx * 3 * (o >= 0 ? 1 : -1),
                  B.y + ny * 3 * (o >= 0 ? 1 : -1));
    cairo_line_to(cr, B.x + nx * ext, B.y + ny * ext);
    cairo_stroke(cr);
    cairo_restore(cr);
    double tw = d->text[0] ? text_width(cr, st, d->text) + 8 : 0;
    double mx = 0.5 * (A2.x + B2.x), my = 0.5 * (A2.y + B2.y);
    if (tw > 0 && tw < L - 24) {
      pil_draw_arrow(cr, pil_v2(mx - dx * tw * 0.5, my - dy * tw * 0.5), A2, 9,
                     st);
      pil_draw_arrow(cr, pil_v2(mx + dx * tw * 0.5, my + dy * tw * 0.5), B2, 9,
                     st);
    } else {
      pil_draw_arrow(cr, pil_v2(mx, my), A2, 9, st);
      pil_draw_arrow(cr, pil_v2(mx, my), B2, 9, st);
    }
    if (tw > 0) {
      pil_draw_text(cr, pil_v2(mx, my), d->text, 0.5, st);
    }
    break;
  }
  case D_ANGLE: {
    double th = d->body ? d->body->angle : 0, r = d->radius;
    cairo_save(cr);
    ink(cr, st);
    cairo_set_line_width(cr, st->thin_w);
    cairo_new_path(cr);
    if (th >= d->ref) {
      cairo_arc(cr, A.x, A.y, r, -th, -d->ref);
    } else {
      cairo_arc(cr, A.x, A.y, r, -d->ref, -th);
    }
    cairo_stroke(cr);
    cairo_restore(cr);
    double mid = -0.5 * (th + d->ref);
    pil_draw_text(cr,
                  pil_v2(A.x + (r + 11) * cos(mid), A.y + (r + 11) * sin(mid)),
                  d->text, 0.5, st);
    break;
  }
  }
}

/* frame */
void pil_draw_world(cairo_t *cr, const pil_world *w, const pil_view *v,
                    const pil_style *st) {
  cairo_save(cr);
  if (st->draw_background) {
    paper(cr, st);
    cairo_paint(cr);
  }

  cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
  cairo_set_antialias(cr, CAIRO_ANTIALIAS_BEST);

  for (int i = 0; i < w->nterr; i++) {
    draw_terrain(cr, w->terrains[i], w->t, v, st);
  }
  for (int i = 0; i < w->ndecos; i++) {
    if (w->decos[i]->kind == D_GROUND) {
      draw_deco(cr, w->decos[i], v, st);
    }
  }
  for (int i = 0; i < w->nels; i++) {
    draw_element(cr, w->els[i], v, st);
  }
  for (int i = 0; i < w->nbodies; i++) {
    draw_body(cr, w->bodies[i], v, st);
  }
  for (int i = 0; i < w->nels; i++) {
    draw_joint_marks(cr, w->els[i], v, st);
  }
  if (st->show_com) {
    for (int i = 0; i < w->nbodies; i++) {
      const pil_body *b = w->bodies[i];
      if (b->shape.kind == pil_SHAPE_POINT ||
          (b->flags & (pil_HIDE_COM | pil_HIDDEN)) || b->dof < 0) {
        continue;
      }
      pil_draw_com(cr, S(v, b->pos), st);
    }
  }
  for (int i = 0; i < w->ndecos; i++)
    if (w->decos[i]->kind != D_GROUND) {
      draw_deco(cr, w->decos[i], v, st);
    }
  for (int i = 0; i < w->nbodies; i++) {
    const pil_body *b = w->bodies[i];
    if (!b->label[0] || (b->flags & pil_HIDDEN)) {
      continue;
    }

    pil_vec2 c = S(v, b->pos);
    double dx = b->label_dx, dy = b->label_dy;
    if (dx == 0 && dy == 0) {
      dx = 14;
      dy = 14;
    }

    pil_draw_text(cr, pil_v2(c.x + dx, c.y - dy), b->label, align_for(dx), st);
  }

  for (int i = 0; i < w->noverlays; i++) {
    w->overlays[i].fn(cr, v, w->overlays[i].user);
  }
  cairo_restore(cr);
}

/* export */
int pil_export_png(const pil_world *w, const pil_view *v, const pil_style *st,
                   const char *path) {
  cairo_surface_t *s = cairo_image_surface_create(
      CAIRO_FORMAT_ARGB32, (int)v->width, (int)v->height);
  cairo_t *cr = cairo_create(s);
  pil_draw_world(cr, w, v, st);
  cairo_status_t r = cairo_surface_write_to_png(s, path);
  cairo_destroy(cr);
  cairo_surface_destroy(s);

  return r == CAIRO_STATUS_SUCCESS ? 0 : -1;
}

#ifdef CAIRO_HAS_PDF_SURFACE
#include <cairo-pdf.h>
#endif
#ifdef CAIRO_HAS_SVG_SURFACE
#include <cairo-svg.h>
#endif

int pil_export_pdf(const pil_world *w, const pil_view *v, const pil_style *st,
                   const char *path) {
#ifdef CAIRO_HAS_PDF_SURFACE
  cairo_surface_t *s = cairo_pdf_surface_create(path, v->width, v->height);
  cairo_t *cr = cairo_create(s);
  pil_draw_world(cr, w, v, st);
  cairo_show_page(cr);
  cairo_destroy(cr);
  cairo_surface_finish(s);
  int ok = cairo_surface_status(s) == CAIRO_STATUS_SUCCESS;
  cairo_surface_destroy(s);
  return ok ? 0 : -1;
#else
  (void)w;
  (void)v;
  (void)st;
  (void)path;
  return -1;
#endif
}

int pil_export_svg(const pil_world *w, const pil_view *v, const pil_style *st,
                   const char *path) {
#ifdef CAIRO_HAS_SVG_SURFACE
  cairo_surface_t *s = cairo_svg_surface_create(path, v->width, v->height);
  cairo_t *cr = cairo_create(s);
  pil_draw_world(cr, w, v, st);
  cairo_destroy(cr);
  cairo_surface_finish(s);
  int ok = cairo_surface_status(s) == CAIRO_STATUS_SUCCESS;
  cairo_surface_destroy(s);
  return ok ? 0 : -1;
#else
  (void)w;
  (void)v;
  (void)st;
  (void)path;
  return -1;
#endif
}

int pil_record(pil_world *w, const pil_view *v, const pil_style *st,
               const char *path, double seconds, double fps, double max_dt) {
  int W = (int)v->width & ~1, H = (int)v->height & ~1;
  size_t plen = strlen(path);
  const char *opts = "-c:v libx264 -pix_fmt yuv420p -crf 18 -preset veryfast";
  if (plen > 4 && !strcmp(path + plen - 4, ".gif")) {
    opts = "-vf "
           "\"split[s0][s1];[s0]palettegen=stats_mode=diff[p];[s1][p]"
           "paletteuse=dither=none\"";
  } else if (plen > 5 && !strcmp(path + plen - 5, ".webm")) {
    opts = "-c:v libvpx-vp9 -pix_fmt yuv420p -crf 28 -b:v 0";
  }

  char cmd[1024];
  snprintf(cmd, sizeof cmd,
           "ffmpeg -loglevel error -y -f rawvideo -pix_fmt bgra -s %dx%d -r %g "
           "-i - %s \"%s\"",
           W, H, fps, opts, path);
  FILE *p = popen(cmd, "w");
  if (!p) {
    return -1;
  }

  cairo_surface_t *s = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, W, H);
  cairo_t *cr = cairo_create(s);
  int frames = (int)(seconds * fps + 0.5);
  for (int f = 0; f < frames; f++) {
    pil_draw_world(cr, w, v, st);
    cairo_surface_flush(s);
    const unsigned char *d = cairo_image_surface_get_data(s);
    int stride = cairo_image_surface_get_stride(s);
    for (int y = 0; y < H; y++) {
      fwrite(d + (size_t)y * stride, 4, (size_t)W, p);
    }

    pil_world_advance((pil_world *)w, 1.0 / fps, max_dt);
  }

  cairo_destroy(cr);
  cairo_surface_destroy(s);

  return pclose(p) == 0 ? 0 : -1;
}
