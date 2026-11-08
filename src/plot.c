/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* plot.c - time-series panels drawn inside figure
 *
 * A plot owns a time axis and up to PIL_PLOT_MAX_SIGNALS curves
 * Curves are computed up front by pil_plots_record() (run world, sample,
 * restore) so panel shows whole result and a cursor follows live animation
 */

#include "plot.h"
#include <cairo.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

pil_plot *pil_add_plot(pil_world *w, double x0, double y0, double x1, double y1,
                       const char *xlabel, const char *ylabel, double t_span) {
  pil_plot *p = calloc(1, sizeof *p);
  if (!p) {
    abort();
  }

  p->x0 = x0;
  p->y0 = y0;
  p->x1 = x1;
  p->y1 = y1;
  p->span = t_span > 0 ? t_span : 1.0;
  snprintf(p->xlabel, sizeof p->xlabel, "%s", xlabel ? xlabel : "");
  snprintf(p->ylabel, sizeof p->ylabel, "%s", ylabel ? ylabel : "");
  if (w->nplots + 1 > w->plcap) {
    w->plcap = (w->nplots + 1) * 2;
    w->plots = realloc(w->plots, (size_t)w->plcap * sizeof *w->plots);
    if (!w->plots) {
      abort();
    }
  }

  w->plots[w->nplots++] = p;

  return p;
}

void pil_plot_add_signal(pil_plot *p, const char *label, pil_signal_fn fn,
                         void *user) {
  if (p->nsig >= PIL_PLOT_MAX_SIGNALS) {
    return;
  }

  pil_signal *s = &p->sig[p->nsig++];
  snprintf(s->label, sizeof s->label, "%s", label ? label : "");
  s->fn = fn;
  s->user = user;
  s->v = NULL;
}

int pil_plot_samples(const pil_plot *p) { return p->n; }

double pil_plot_time(const pil_plot *p, int k) { return p->t[k]; }

double pil_plot_value(const pil_plot *p, int sig, int k) {
  return p->sig[sig].v[k];
}

void pil_plots_free(pil_world *w) {
  for (int i = 0; i < w->nplots; i++) {
    pil_plot *p = w->plots[i];
    for (int j = 0; j < p->nsig; j++) {
      free(p->sig[j].v);
    }

    free(p->t);
    free(p);
  }

  free(w->plots);
  w->plots = NULL;
  w->nplots = 0;
}

void pil_plots_record(pil_world *w, double sample_dt, double max_dt) {
  if (w->nplots == 0 || sample_dt <= 0) {
    return;
  }

  double span = 0;
  for (int i = 0; i < w->nplots; i++) {
    span = fmax(span, w->plots[i]->span);
  }

  int n = (int)floor(span / sample_dt + 1e-9) + 1;

  // Remember state: run must leave figure exactly as it was
  int nb = w->nbodies;
  double *save = malloc((size_t)nb * 6 * sizeof *save);
  if (!save) {
    abort();
  }

  for (int i = 0; i < nb; i++) {
    const pil_body *b = w->bodies[i];

    save[6 * i + 0] = b->pos.x;
    save[6 * i + 1] = b->pos.y;
    save[6 * i + 2] = b->angle;
    save[6 * i + 3] = b->vel.x;
    save[6 * i + 4] = b->vel.y;
    save[6 * i + 5] = b->omega;
  }

  double t0 = w->t;
  double diss0 = w->dissipated;

  for (int i = 0; i < w->nplots; i++) {
    pil_plot *p = w->plots[i];
    free(p->t);
    p->t = malloc((size_t)n * sizeof *p->t);
    for (int j = 0; j < p->nsig; j++) {
      free(p->sig[j].v);
      p->sig[j].v = malloc((size_t)n * sizeof *p->sig[j].v);
    }

    p->n = n;
    p->cap = n;
  }

  for (int k = 0; k < n; k++) {
    for (int i = 0; i < w->nplots; i++) {
      pil_plot *p = w->plots[i];
      p->t[k] = k * sample_dt;
      for (int j = 0; j < p->nsig; j++) {
        p->sig[j].v[k] = p->sig[j].fn(w, p->sig[j].user);
      }
    }

    if (k + 1 < n) {
      pil_world_advance(w, sample_dt, max_dt);
    }
  }

  for (int i = 0; i < nb; i++) {
    pil_body *b = w->bodies[i];
    b->pos = pil_v2(save[6 * i + 0], save[6 * i + 1]);
    b->angle = save[6 * i + 2];
    b->vel = pil_v2(save[6 * i + 3], save[6 * i + 4]);
    b->omega = save[6 * i + 5];
  }

  free(save);
  w->t = t0;
  w->tcur = t0;
  w->dissipated = diss0;
  for (int i = 0; i < w->nels; i++) {
    w->els[i]->rope_off = 0;
  }

  pil_elements_rebase(w);
  pil_circuits_reset(w); // recording starts at t = 0
  pil_lsys_reset_all(w);
}

// Drawing

static double nice_step(double range, int target) {
  double raw = range / target;
  double mag = pow(10.0, floor(log10(raw)));
  double n = raw / mag;

  return mag * (n < 1.5 ? 1.0 : (n < 3.0 ? 2.0 : (n < 7.0 ? 5.0 : 10.0)));
}

static void ink(cairo_t *cr, const pil_style *st) {
  cairo_set_source_rgb(cr, st->ink[0], st->ink[1], st->ink[2]);
}

static void number(char *buf, size_t n, double v, double step) {
  if (fabs(v) < step * 1e-6) {
    v = 0;
  }

  snprintf(buf, n, "%g", v);
}

static double sample_at(const pil_plot *p, int j, double t) {
  if (p->n < 2) {
    return p->n ? p->sig[j].v[0] : 0.0;
  }

  double dt = p->t[1] - p->t[0], f = t / dt;
  int k = (int)floor(f);
  if (k < 0) {
    return p->sig[j].v[0];
  }

  if (k >= p->n - 1) {
    return p->sig[j].v[p->n - 1];
  }

  f -= k;

  return p->sig[j].v[k] * (1 - f) + p->sig[j].v[k + 1] * f;
}

void pil_plot_draw(cairo_t *cr, const pil_plot *p, const pil_world *w,
                   const pil_view *v, const pil_style *st) {
  const double ml = 66;
  const double mr = 14;
  const double mt = 10;
  const double mb = 38;
  double bx0 = p->x0 * v->width + ml, bx1 = p->x1 * v->width - mr;
  double by0 = p->y0 * v->height + mt, by1 = p->y1 * v->height - mb;
  if (bx1 - bx0 < 40 || by1 - by0 < 30) {
    return;
  }

  // y range from data
  double lo = 0;
  double hi = 1;
  int have = 0;
  for (int j = 0; j < p->nsig && p->n > 0; j++) {
    for (int k = 0; k < p->n; k++) {
      double y = p->sig[j].v[k];
      if (!have) {
        lo = hi = y;
        have = 1;
      }

      lo = fmin(lo, y);
      hi = fmax(hi, y);
    }
  }

  if (hi - lo < 1e-12) {
    hi = lo + 1.0;
  }

  // Headroom for legend, then snap to tick marks (a value within 3% of a tick
  // sits on it instead of forcing a whole extra step)
  double legend_frac = (p->nsig * 16 + 10) / (by1 - by0);
  hi += (hi - lo) * legend_frac / (1.0 - legend_frac);
  double ystep = nice_step(hi - lo, 4), snap = 0.03 * (hi - lo);
  lo = floor((lo + snap) / ystep) * ystep;
  hi = ceil((hi - snap) / ystep) * ystep;
  double xstep = nice_step(p->span, 5);

#define PX(t) (bx0 + (t) / p->span * (bx1 - bx0))
#define PY(y) (by1 - ((y) - lo) / (hi - lo) * (by1 - by0))

  pil_style up = *st; /* digits stand upright, symbols stay italic */
  up.upright = 1;
  cairo_save(cr);
  ink(cr, st);
  cairo_set_line_width(cr, st->thin_w);
  cairo_rectangle(cr, bx0, by0, bx1 - bx0, by1 - by0);
  cairo_stroke(cr);

  char buf[32];
  for (double t = 0; t <= p->span + 1e-9; t += xstep) {
    cairo_move_to(cr, PX(t), by1);
    cairo_line_to(cr, PX(t), by1 - 4);
    cairo_stroke(cr);
    number(buf, sizeof buf, t, xstep);
    pil_draw_text(cr, pil_v2(PX(t), by1 + 14), buf, 0.5, &up);
  }

  for (double y = lo; y <= hi + ystep * 1e-6; y += ystep) {
    cairo_move_to(cr, bx0, PY(y));
    cairo_line_to(cr, bx0 + 4, PY(y));
    cairo_stroke(cr);
    number(buf, sizeof buf, y, ystep);
    pil_draw_text(cr, pil_v2(bx0 - 6, PY(y) + 4), buf, 1.0, &up);
  }

  pil_draw_text(cr, pil_v2((bx0 + bx1) / 2, by1 + 31), p->xlabel, 0.5, st);
  cairo_save(cr);
  cairo_translate(cr, p->x0 * v->width + 12, (by0 + by1) / 2);
  cairo_rotate(cr, -M_PI / 2);
  pil_draw_text(cr, pil_v2(0, 4), p->ylabel, 0.5, st);
  cairo_restore(cr);

  // Curves: distinct dash patterns, all in ink
  static const double dashes[5][4] = {{0, 0, 0, 0},
                                      {7, 3, 7, 3},
                                      {1.5, 3, 1.5, 3},
                                      {10, 3, 2, 3},
                                      {4, 2, 1, 2}};
  cairo_save(cr);
  cairo_rectangle(cr, bx0, by0, bx1 - bx0, by1 - by0);
  cairo_clip(cr);
  for (int j = 0; j < p->nsig && p->n > 1; j++) {
    int style = j % 5;
    cairo_set_line_width(cr, st->thin_w * 1.15);
    if (style == 0) {
      cairo_set_dash(cr, NULL, 0, 0);
    } else {
      cairo_set_dash(cr, dashes[style], 4, 0);
    }

    for (int k = 0; k < p->n; k++) {
      if (k == 0) {
        cairo_move_to(cr, PX(p->t[k]), PY(p->sig[j].v[k]));
      } else {
        cairo_line_to(cr, PX(p->t[k]), PY(p->sig[j].v[k]));
      }
    }

    cairo_stroke(cr);
  }

  cairo_restore(cr);

  // Cursor and markers at current time
  if (p->n > 1 && w->t >= 0 && w->t <= p->span) {
    cairo_set_dash(cr, NULL, 0, 0);
    cairo_set_line_width(cr, st->thin_w * 0.7);
    cairo_move_to(cr, PX(w->t), by0);
    cairo_line_to(cr, PX(w->t), by1);
    cairo_stroke(cr);
    for (int j = 0; j < p->nsig; j++) {
      cairo_arc(cr, PX(w->t), PY(sample_at(p, j, w->t)), 3.0, 0, 2 * M_PI);
      cairo_fill(cr);
    }
  }

  // Legend, top right inside box
  double ly = by0 + 14;
  for (int j = 0; j < p->nsig; j++) {
    if (!p->sig[j].label[0]) {
      continue;
    }

    int style = j % 5;
    cairo_set_line_width(cr, st->thin_w * 1.15);
    if (style == 0) {
      cairo_set_dash(cr, NULL, 0, 0);
    } else {
      cairo_set_dash(cr, dashes[style], 4, 0);
    }

    cairo_move_to(cr, bx1 - 40, ly - 4);
    cairo_line_to(cr, bx1 - 12, ly - 4);
    cairo_stroke(cr);
    cairo_set_dash(cr, NULL, 0, 0);
    pil_draw_text(cr, pil_v2(bx1 - 48, ly), p->sig[j].label, 1.0, st);
    ly += 16;
  }

  cairo_restore(cr);
#undef PX
#undef PY
}
