/*
 * hohmann.c - Hohmann transfer between circular orbits, r1 = 1 to r2 = 2 (mu =
 * 1) Two tangential burns: dv1 = sqrt(mu/r1) (sqrt(2 r2 / (r1 + r2)) - 1) puts
 * craft on transfer ellipse with periapsis r1 and apoapsis r2; half an orbit
 * later, dv2 = sqrt(mu/r2) (1 - sqrt(2 r1 / (r1 + r2))) circularises it
 * Flight time pi sqrt(a^3 / mu) with a = (r1 + r2) / 2
 * Plots: radius and specific energy -mu / 2a against time, with energy steps of
 * burns
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 11.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_orbit.h>

#define R1 1.0
#define R2 2.0
#define DT 0.0005
#define TEND 20.0
#define T1 1.0

typedef struct {
  double *x;
  double *y;
  double *r;
  double *e;
  int n;
  double tof;
  const pil_world *w;
} hoh_t;

static int idx(const hoh_t *h) {
  int i = (int)(pil_world_time(h->w) / DT);

  return i < h->n ? i : h->n - 1;
}

static void circle(cairo_t *cr, pil_vec2 c, double r, int dashed) {
  double dash[2] = {6, 5};
  if (dashed) {
    cairo_set_dash(cr, dash, 2, 0);
  }
  cairo_arc(cr, c.x, c.y, r, 0, 2 * M_PI);
  cairo_stroke(cr);
  cairo_set_dash(cr, NULL, 0, 0);
}

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const hoh_t *h = user;
  pil_style st = pil_style_textbook();
  int cur = idx(h);
  pil_vec2 O = pil_view_to_screen(v, pil_v2(0, 0));
  double sc = pil_view_to_screen(v, pil_v2(1, 0)).x - O.x;

  cairo_new_path(cr);
  cairo_set_line_width(cr, st.thin_w);
  cairo_set_source_rgb(cr, 0.55, 0.55, 0.55);
  circle(cr, O, R1 * sc, 0);
  circle(cr, O, R2 * sc, 0);

  // Transfer ellipse (dashed): centre at -(a - r1) on x axis, semi-axes a, b
  double a = 0.5 * (R1 + R2);
  double e = (R2 - R1) / (R1 + R2);
  double b = a * sqrt(1 - e * e);
  cairo_save(cr);
  cairo_translate(cr, O.x, O.y);
  cairo_rotate(
      cr, -T1 * sqrt(1.0 / (R1 * R1 * R1))); // Burn 1 happens T1 after start
  cairo_translate(cr, -(a - R1) * sc, 0);
  cairo_scale(cr, a * sc, b * sc);
  cairo_arc(cr, 0, 0, 1.0, 0, 2 * M_PI);
  cairo_restore(cr);
  double dash[2] = {6, 5};
  cairo_set_dash(cr, dash, 2, 0);
  cairo_set_line_width(cr, st.thin_w);
  cairo_stroke(cr);
  cairo_set_dash(cr, NULL, 0, 0);

  // Trail so far
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_set_line_width(cr, st.line_w);
  for (int i = 0; i <= cur; i += 20) {
    pil_vec2 q = pil_view_to_screen(v, pil_v2(h->x[i], h->y[i]));
    if (i == 0) {
      cairo_move_to(cr, q.x, q.y);
    } else {
      cairo_line_to(cr, q.x, q.y);
    }
  }

  cairo_stroke(cr);
  pil_vec2 q = pil_view_to_screen(v, pil_v2(h->x[cur], h->y[cur]));
  cairo_arc(cr, q.x, q.y, 6, 0, 2 * M_PI);
  cairo_fill(cr);

  // Central body and burn markers
  cairo_arc(cr, O.x, O.y, 9, 0, 2 * M_PI);
  cairo_fill(cr);
  double ph = T1 * sqrt(1.0 / (R1 * R1 * R1));
  pil_vec2 b1 =
      pil_view_to_screen(v, pil_v2(1.12 * R1 * cos(ph), 1.12 * R1 * sin(ph)));
  pil_vec2 b2 =
      pil_view_to_screen(v, pil_v2(-1.12 * R2 * cos(ph), -1.12 * R2 * sin(ph)));
  pil_draw_text(cr, b1, "burn 1", 0.5, &st);
  pil_draw_text(cr, b2, "burn 2", 0.5, &st);
  cairo_new_path(cr);
}

static double sig_r(const pil_world *w, void *u) {
  (void)w;
  const hoh_t *h = u;

  return h->r[idx(h)];
}

static double sig_e(const pil_world *w, void *u) {
  (void)w;
  const hoh_t *h = u;

  return h->e[idx(h)];
}

int main(int argc, char **argv) {
  static hoh_t h;
  h.n = (int)(TEND / DT);
  h.x = malloc((size_t)h.n * sizeof *h.x);
  h.y = malloc((size_t)h.n * sizeof *h.y);
  h.r = malloc((size_t)h.n * sizeof *h.r);
  h.e = malloc((size_t)h.n * sizeof *h.e);
  double at = 0.5 * (R1 + R2);
  h.tof = M_PI * sqrt(at * at * at);
  double dv1 = sqrt(1.0 / R1) * (sqrt(2 * R2 / (R1 + R2)) - 1);
  double dv2 = sqrt(1.0 / R2) * (1 - sqrt(2 * R1 / (R1 + R2)));

  // Start on x axis in a circular orbit, counter-clockwise
  double x = R1;
  double y = 0;
  double vx = 0;
  double vy = sqrt(1.0 / R1);
  int burn1 = (int)(T1 / DT + 0.5);
  double burn2 = (int)((T1 + h.tof) / DT + 0.5);
  for (int i = 0; i < h.n; i++) {
    double r = hypot(x, y);
    h.x[i] = x;
    h.y[i] = y;
    h.r[i] = r;
    h.e[i] = 0.5 * (vx * vx + vy * vy) - 1.0 / r;
    if (i == burn1 || i == burn2) {
      double sp = hypot(vx, vy);
      double dv = i == burn1 ? dv1 : dv2;
      vx += dv * vx / sp;
      vy += dv * vy / sp;
    }

    // Kick-drift-kick leapfrog
    double r3 = r * r * r;
    vx -= 0.5 * DT * x / r3;
    vy -= 0.5 * DT * y / r3;
    x += DT * vx;
    y += DT * vy;
    r = hypot(x, y);
    r3 = r * r * r;
    vx -= 0.5 * DT * x / r3;
    vy -= 0.5 * DT * y / r3;
  }

  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  h.w = w;
  pil_world_add_overlay(w, overlay, &h);

  pil_plot *p = pil_add_plot(w, 0.62, 0.06, 0.99, 0.46, "t", "r", TEND);
  pil_plot_add_signal(p, "r", sig_r, &h);
  pil_plot *q = pil_add_plot(w, 0.62, 0.56, 0.99, 0.97, "t", "energy", TEND);
  pil_plot_add_signal(q, "E", sig_e, &h);
  pil_plots_record(w, 0.05, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -2.6, -2.2, 2.6, 2.2, 1300, 760, 0.0, 0.0, 0.6,
                      1.0, 14);
  int rc = example_main(w, &view, "Hohmann transfer orbit", argc, argv,
                        "hohmann.png");

  pil_world_free(w);
  free(h.x);
  free(h.y);
  free(h.r);
  free(h.e);

  return rc;
}
