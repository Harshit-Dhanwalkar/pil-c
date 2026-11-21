/*
 * foucault.c - Foucault pendulum seen from above, in rotating frame of Earth
 * In frame rotating at angular rate Omega pendulum feels
 * Coriolis force, so for small swings
 *     x'' = -w0^2 x + 2 W y',   y'' = -w0^2 y - 2 W x',   W = Omega
 * sin(latitude) and swing plane turns at rate -W (clockwise in northern
 * hemisphere), a full turn in one sidereal day / sin(latitude)
 * W is exaggerated here so precession is visible in seconds
 * Right: plane angle measured from path over each swing, against prediction
 * -W t
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 30.0
#include "example_main.h"
#include <math.h>

#define W0 6.283185307179586
#define WR 0.15
#define DT 0.002
#define TEND 60.0

typedef struct {
  double *x;
  double *y;
  double *psi;
  int n;
  const pil_world *w;
} fou_t;

static void deriv(const double s[4], double d[4]) {
  d[0] = s[2];
  d[1] = s[3];
  d[2] = -W0 * W0 * s[0] + 2 * WR * s[3];
  d[3] = -W0 * W0 * s[1] - 2 * WR * s[2];
}

static void rk4(double s[4], double h) {
  double k1[4];
  double k2[4];
  double k3[4];
  double k4[4];
  double t[4];
  deriv(s, k1);
  for (int i = 0; i < 4; i++) {
    t[i] = s[i] + 0.5 * h * k1[i];
  }
  deriv(t, k2);
  for (int i = 0; i < 4; i++) {
    t[i] = s[i] + 0.5 * h * k2[i];
  }
  deriv(t, k3);
  for (int i = 0; i < 4; i++) {
    t[i] = s[i] + h * k3[i];
  }
  deriv(t, k4);
  for (int i = 0; i < 4; i++) {
    s[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
  }
}

static int idx(const fou_t *f) {
  int i = (int)(pil_world_time(f->w) / DT);

  return i < f->n ? i : f->n - 1;
}

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const fou_t *f = user;
  pil_style st = pil_style_textbook();
  int cur = idx(f);
  pil_vec2 O = pil_view_to_screen(v, pil_v2(0, 0));

  // Axes
  cairo_new_path(cr);
  cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
  cairo_set_line_width(cr, st.thin_w);
  pil_vec2 a = pil_view_to_screen(v, pil_v2(-1.3, 0));
  pil_vec2 b = pil_view_to_screen(v, pil_v2(1.3, 0));
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, b.x, b.y);
  a = pil_view_to_screen(v, pil_v2(0, -1.3));
  b = pil_view_to_screen(v, pil_v2(0, 1.3));
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, b.x, b.y);
  cairo_stroke(cr);

  // Predicted swing plane (dashed)
  double th = -WR * pil_world_time(f->w);
  double dash[2] = {6, 5};
  cairo_set_dash(cr, dash, 2, 0);
  a = pil_view_to_screen(v, pil_v2(-1.2 * cos(th), -1.2 * sin(th)));
  b = pil_view_to_screen(v, pil_v2(1.2 * cos(th), 1.2 * sin(th)));
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, b.x, b.y);
  cairo_stroke(cr);
  cairo_set_dash(cr, NULL, 0, 0);

  // Wake: last 8 s, fading
  int span = (int)(8.0 / DT);
  int first = cur - span > 0 ? cur - span : 0;
  cairo_set_line_width(cr, 1.2);
  for (int i = first + 1; i <= cur; i++) {
    double g = 0.05 + 0.9 * (double)(cur - i) / span;
    pil_vec2 p0 = pil_view_to_screen(v, pil_v2(f->x[i - 1], f->y[i - 1]));
    pil_vec2 p1 = pil_view_to_screen(v, pil_v2(f->x[i], f->y[i]));
    cairo_set_source_rgb(cr, g, g, g);
    cairo_move_to(cr, p0.x, p0.y);
    cairo_line_to(cr, p1.x, p1.y);
    cairo_stroke(cr);
  }

  pil_vec2 q = pil_view_to_screen(v, pil_v2(f->x[cur], f->y[cur]));
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_arc(cr, q.x, q.y, 7, 0, 2 * M_PI);
  cairo_fill(cr);
  cairo_arc(cr, O.x, O.y, 3, 0, 2 * M_PI);
  cairo_fill(cr);
  cairo_new_path(cr);
  (void)st;
}

static double sig_meas(const pil_world *w, void *u) {
  (void)w;
  const fou_t *f = u;

  return f->psi[idx(f)];
}

static double sig_theory(const pil_world *w, void *u) {
  (void)u;

  return -WR * pil_world_time(w);
}

int main(int argc, char **argv) {
  static fou_t f;
  f.n = (int)(TEND / DT);
  f.x = malloc((size_t)f.n * sizeof *f.x);
  f.y = malloc((size_t)f.n * sizeof *f.y);
  f.psi = malloc((size_t)f.n * sizeof *f.psi);

  double s[4] = {1.0, 0, 0, 0}; // Released from rest at x = 1
  for (int i = 0; i < f.n; i++) {
    f.x[i] = s[0];
    f.y[i] = s[1];
    rk4(s, DT);
  }

  // Plane angle from second moments of last full period, unwrapped
  int per = (int)(1.0 / DT);
  double prev = 0;
  for (int i = 0; i < f.n; i++) {
    int j0 = i - per > 0 ? i - per : 0;
    double sxx = 0;
    double syy = 0;
    double sxy = 0;
    for (int j = j0; j <= i; j++) {
      sxx += f.x[j] * f.x[j];
      syy += f.y[j] * f.y[j];
      sxy += f.x[j] * f.y[j];
    }

    double ps = i < per / 2 ? 0.0 : 0.5 * atan2(2 * sxy, sxx - syy);
    while (ps - prev > M_PI / 2) {
      ps -= M_PI;
    }
    while (ps - prev < -M_PI / 2) {
      ps += M_PI;
    }
    f.psi[i] = ps;
    prev = ps;
  }

  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  f.w = w;
  pil_world_add_overlay(w, overlay, &f);

  pil_plot *p = pil_add_plot(w, 0.62, 0.20, 0.99, 0.80, "t [s]",
                             "plane angle [rad]", 30.0);
  pil_plot_add_signal(p, "measured", sig_meas, &f);
  pil_plot_add_signal(p, "-W t", sig_theory, &f);
  pil_plots_record(w, 0.1, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -1.4, -1.4, 1.4, 1.4, 1300, 760, 0.0, 0.0, 0.6,
                      1.0, 14);
  int rc = example_main(w, &view, "Foucault pendulum (top view)", argc, argv,
                        "foucault.png");

  pil_world_free(w);
  free(f.x);
  free(f.y);
  free(f.psi);

  return rc;
}
