/*
 * poincare.c - driven damped pendulum and its Poincare section
 *     theta'' + theta'/q + sin(theta) = g cos(w_D t)
 *     q = 2, w_D = 2/3, g = 1.5
 * a route to chaos
 * Grey curve is phase-space trajectory (theta wrapped to -pi..pi, omega); black
 * dots are stroboscopic samples taken once per drive period
 * A regular motion would leave a few dots; here they build up a fractal strange
 * attractor
 * Time runs fast: one second of animation is 40 drive-time units
 * Everything is public-API
 */

#include <stddef.h>
#include <stdlib.h>
#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 250.0
#include "example_main.h"
#include <math.h>

#define Q 2.0
#define WD (2.0 / 3.0)
#define G 1.5
#define DT 0.01
#define RATE 40.0
#define NPER 4000

typedef struct {
  double *th;      // Wrapped angle at every step
  double *om;      // Wrapped velocity at every step
  double (*sp)[2]; // Strobe points
  int n;
  int ns;
  const pil_world *w;
} poin_t;

static void deriv(double t, const double s[2], double d[2]) {
  d[0] = s[1];
  d[1] = -s[1] / Q - sin(s[0]) + G * cos(WD * t);
}

static void rk4(double t, double s[2], double h) {
  double k1[2];
  double k2[2];
  double k3[2];
  double k4[2];
  double u[2];
  deriv(t, s, k1);
  u[0] = s[0] + 0.5 * h * k1[0];
  u[1] = s[1] + 0.5 * h * k1[1];
  deriv(t + 0.5 * h, u, k2);
  u[0] = s[0] + 0.5 * h * k2[0];
  u[1] = s[1] + 0.5 * h * k2[1];
  deriv(t + 0.5 * h, u, k3);
  u[0] = s[0] + h * k3[0];
  u[1] = s[1] + h * k3[1];
  deriv(t + h, u, k4);
  for (int i = 0; i < 2; i++) {
    s[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
}
}

static double wrap(double a) {
  a = fmod(a + M_PI, 2 * M_PI);
  if (a < 0) {
    a += 2 * M_PI;
}

  return a - M_PI;
}

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const poin_t *p = user;
  pil_style st = pil_style_textbook();
  int cur = (int)(pil_world_time(p->w) * RATE / DT);
  if (cur >= p->n) {
    cur = p->n - 1;
}

  // Axes
  cairo_new_path(cr);
  cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
  cairo_set_line_width(cr, st.thin_w);
  pil_vec2 a = pil_view_to_screen(v, pil_v2(-M_PI, 0)),
           b = pil_view_to_screen(v, pil_v2(M_PI, 0));
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, b.x, b.y);
  a = pil_view_to_screen(v, pil_v2(0, -2.6));
  b = pil_view_to_screen(v, pil_v2(0, 2.6));
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, b.x, b.y);
  cairo_stroke(cr);
  cairo_new_path(cr);
  cairo_rectangle(cr, pil_view_to_screen(v, pil_v2(-M_PI, 2.6)).x,
                  pil_view_to_screen(v, pil_v2(-M_PI, 2.6)).y,
                  pil_view_to_screen(v, pil_v2(M_PI, 0)).x -
                      pil_view_to_screen(v, pil_v2(-M_PI, 0)).x,
                  pil_view_to_screen(v, pil_v2(0, -2.6)).y -
                      pil_view_to_screen(v, pil_v2(0, 2.6)).y);
  cairo_stroke(cr);
  cairo_new_path(cr);

  // Trajectory of the last 25 drive-time units
  int first = cur - (int)(25.0 / DT) > 0 ? cur - (int)(25.0 / DT) : 0;
  cairo_set_line_width(cr, 0.8);
  for (int i = first + 1; i <= cur; i++) {
    if (fabs(p->th[i] - p->th[i - 1]) > 3.0) {
      continue; // Wrap jump
      }
    double age = (double)(cur - i) / (cur - first + 1); double g = 0.2 + 0.7 * age;
    pil_vec2 p0 = pil_view_to_screen(v, pil_v2(p->th[i - 1], p->om[i - 1]));
    pil_vec2 p1 = pil_view_to_screen(v, pil_v2(p->th[i], p->om[i]));
    cairo_set_source_rgb(cr, g, g, g);
    cairo_move_to(cr, p0.x, p0.y);
    cairo_line_to(cr, p1.x, p1.y);
    cairo_stroke(cr);
  }

  // Strobe points revealed so far
  double T = 2 * M_PI / WD;
  int ns = (int)(cur * DT / T);
  if (ns > p->ns) {
    ns = p->ns;}
  cairo_set_source_rgb(cr, 0, 0, 0);
  for (int i = 0; i < ns; i++) {
    pil_vec2 q = pil_view_to_screen(v, pil_v2(p->sp[i][0], p->sp[i][1]));
    cairo_rectangle(cr, q.x - 1.0, q.y - 1.0, 2.0, 2.0);
  }

  cairo_fill(cr);
  pil_vec2 q = pil_view_to_screen(v, pil_v2(p->th[cur], p->om[cur]));
  cairo_arc(cr, q.x, q.y, 5, 0, 2 * M_PI);
  cairo_fill(cr);
  pil_draw_text(cr, pil_view_to_screen(v, pil_v2(0, -2.85)),
                "\\theta  (wrapped)", 0.5, &st);
  pil_draw_text(cr, pil_view_to_screen(v, pil_v2(-M_PI - 0.25, 0)), "\\omega",
                1.0, &st);
  cairo_new_path(cr);
}

int main(int argc, char **argv) {
  static poin_t p;
  double T = 2 * M_PI / WD;
  p.sp = malloc(NPER * sizeof *p.sp);

  double s[2] = {0.2, 0.0};
  int per = (int)lround(T / DT); // Steps per drive period
  double h = T / per;
  p.n = NPER * per;
  p.th = malloc((size_t)p.n * sizeof *p.th);
  p.om = malloc((size_t)p.n * sizeof *p.om);
  for (int i = 0; i < p.n; i++) {
    p.th[i] = wrap(s[0]);
    p.om[i] = s[1];
    if (i % per == 0) {
      p.sp[i / per][0] = wrap(s[0]);
      p.sp[i / per][1] = s[1];
    }

    rk4(i * h, s, h);
  }

  p.ns = NPER;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  p.w = w;
  pil_world_add_overlay(w, overlay, &p);

  pil_view view;
  pil_view_fit_region(&view, -M_PI - 0.6, -3.1, M_PI + 0.6, 3.1, 1300, 760, 0.0,
                      0.0, 1.0, 1.0, 14);
  int rc = example_main(w, &view, "Driven pendulum: Poincare section", argc,
                        argv, "poincare.png");

  pil_world_free(w);
  free(p.th);
  free(p.om);
  free(p.sp);

  return rc;
}
