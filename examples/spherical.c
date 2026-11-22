/*
 * spherical.c - spherical pendulum from its Lagrangian
 * Bob on a rod of length l, polar angle theta from downward vertical and
 * azimuth phi:
 *   T = 1/2 m l^2 (th'^2 + sin^2 th phi'^2)
 *   V = -m g l cos th
 * phi is cyclic, so L_z = m l^2 sin^2 th phi' is conserved (plotted flat)
 * while theta oscillates between two turning angles: bob's shadow on floor
 * traces a rosette
 * Started with a little azimuthal speed, below value 'phi'^2 = g / (l cos th)'
 * of circular (conical) pendulum
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 14.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_lagrange.h>

#define TRACE_DT 0.02

typedef struct {
  double m, l, g;
} sp_params;

typedef struct {
  const sp_params *p;
  pil_lsys *s;
  const pil_world *w;
  pil_proj3 pr;
  double (*trace)[3];
  int ntrace;
} sp_scene;

static double sp_T(const double *q, const double *qd, void *user) {
  const sp_params *p = user;
  double s = sin(q[0]);

  return 0.5 * p->m * p->l * p->l * (qd[0] * qd[0] + s * s * qd[1] * qd[1]);
}

static double sp_V(const double *q, void *user) {
  const sp_params *p = user;

  return -p->m * p->g * p->l * cos(q[0]);
}

// Bob position, z up, pivot at origin
static void bob(const sp_params *p, const double *q, double *x, double *y,
                double *z) {
  *x = p->l * sin(q[0]) * cos(q[1]);
  *y = p->l * sin(q[0]) * sin(q[1]);
  *z = -p->l * cos(q[0]);
}

static void sp_overlay(cairo_t *cr, const pil_view *v, void *user) {
  const sp_scene *sc = user;
  const sp_params *p = sc->p;
  pil_style st = pil_style_textbook();
  double zf = -1.18 * p->l; // floor plane
  double t = pil_world_time(sc->w);

  // Floor square and axes
  cairo_set_source_rgb(cr, 0.55, 0.55, 0.55);
  cairo_set_line_width(cr, st.thin_w);
  double h = 1.1 * p->l;
  pil_vec2 c[4];
  double sx[4] = {-1, 1, 1, -1};
  double sy[4] = {-1, -1, 1, 1};
  for (int i = 0; i < 4; i++) {
    c[i] = pil_view_to_screen(
        v, pil_proj3_point(&sc->pr, sx[i] * h, sy[i] * h, zf));
  }

  cairo_move_to(cr, c[0].x, c[0].y);
  for (int i = 1; i < 4; i++) {
    cairo_line_to(cr, c[i].x, c[i].y);
  }

  cairo_close_path(cr);
  cairo_stroke(cr);

  // Shadow rosette (path so far, fading), and bob's own wake
  int last = (int)(t / TRACE_DT);
  if (last >= sc->ntrace) {
    last = sc->ntrace - 1;
  }

  for (int pass = 0; pass < 2; pass++) {
    for (int i = 1; i <= last; i++) {
      double zz = pass ? sc->trace[i][2] : zf;
      pil_vec2 a = pil_view_to_screen(
          v, pil_proj3_point(&sc->pr, sc->trace[i - 1][0], sc->trace[i - 1][1],
                             pass ? sc->trace[i - 1][2] : zf));
      pil_vec2 b = pil_view_to_screen(
          v, pil_proj3_point(&sc->pr, sc->trace[i][0], sc->trace[i][1], zz));
      double g = pass ? 0.35 : 0.0;
      cairo_set_source_rgb(cr, g, g, g);
      cairo_set_line_width(cr, pass ? 0.9 : 1.1);
      cairo_move_to(cr, a.x, a.y);
      cairo_line_to(cr, b.x, b.y);
      cairo_stroke(cr);
    }
  }

  double x;
  double y;
  double z;
  bob(p, pil_lsys_q(sc->s), &x, &y, &z);
  pil_vec2 O = pil_view_to_screen(v, pil_proj3_point(&sc->pr, 0, 0, 0));
  pil_vec2 B = pil_view_to_screen(v, pil_proj3_point(&sc->pr, x, y, z));
  pil_vec2 S = pil_view_to_screen(v, pil_proj3_point(&sc->pr, x, y, zf));
  pil_vec2 F = pil_view_to_screen(v, pil_proj3_point(&sc->pr, 0, 0, zf));
  pil_vec2 Zt =
      pil_view_to_screen(v, pil_proj3_point(&sc->pr, 0, 0, 0.2 * p->l));

  // Vertical axis (dash-dot), plumb line to shadow, rod and bob
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_set_line_width(cr, st.thin_w);
  double dash[2] = {8, 3};
  cairo_set_dash(cr, dash, 2, 0);
  cairo_move_to(cr, Zt.x, Zt.y);
  cairo_line_to(cr, F.x, F.y);
  cairo_move_to(cr, B.x, B.y);
  cairo_line_to(cr, S.x, S.y);
  cairo_stroke(cr);
  cairo_set_dash(cr, NULL, 0, 0);
  cairo_set_line_width(cr, st.line_w);
  cairo_move_to(cr, O.x, O.y);
  cairo_line_to(cr, B.x, B.y);
  cairo_stroke(cr);
  cairo_arc(cr, B.x, B.y, 9, 0, 2 * M_PI);
  cairo_set_source_rgb(cr, 1, 1, 1);
  cairo_fill_preserve(cr);
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_stroke(cr);
  cairo_arc(cr, S.x, S.y, 3, 0, 2 * M_PI);
  cairo_fill(cr);
  cairo_arc(cr, O.x, O.y, 3.5, 0, 2 * M_PI);
  cairo_fill(cr);
  pil_draw_text(cr, pil_v2(B.x + 14, B.y - 4), "m", 0.0, &st);
  pil_draw_text(cr, pil_v2(O.x + 10, O.y - 10), "O", 0.0, &st);
  pil_draw_text(cr, pil_v2(Zt.x - 6, Zt.y - 8), "z", 0.5, &st);
}

static double sig_theta(const pil_world *w, void *u) {
  (void)w;

  return pil_lsys_q(u)[0];
}

typedef struct {
  const sp_params *p;
  pil_lsys *s;
} lz_ref;

static double sig_lz(const pil_world *w, void *u) {
  (void)w;
  const lz_ref *r = u;
  double s = sin(pil_lsys_q(r->s)[0]);

  return r->p->m * r->p->l * r->p->l * s * s * pil_lsys_qd(r->s)[1];
}

static double sig_E(const pil_world *w, void *u) {
  (void)w;

  return pil_lsys_jacobi(u);
}

int main(int argc, char **argv) {
  static sp_params p = {1.0, 1.0, 9.81};
  static sp_scene sc;
  static lz_ref lz;
  const double q0[2] = {0.9, 0.0};
  const double qd0[2] = {0.0, 2.2};

  {
    pil_world *w2 = pil_world_new();
    pil_lsys *s2 = pil_add_lagrangian(w2, 2, sp_T, sp_V, &p);
    pil_lsys_set_state(s2, q0, qd0);
    int n = (int)(40.0 / TRACE_DT);
    sc.trace = malloc((size_t)n * sizeof *sc.trace);
    sc.ntrace = n;
    for (int i = 0; i < n; i++) {
      bob(&p, pil_lsys_q(s2), &sc.trace[i][0], &sc.trace[i][1],
          &sc.trace[i][2]);
      pil_world_advance(w2, TRACE_DT, 1e-3);
    }

    pil_world_free(w2);
  }

  pil_world *w = pil_world_new();
  sc.p = &p;
  sc.w = w;
  sc.pr = pil_proj3_make(0, 0, 0.5, 0.45, 1.0);
  sc.s = pil_add_lagrangian(w, 2, sp_T, sp_V, &p);
  pil_lsys_set_state(sc.s, q0, qd0);
  lz.p = &p;
  lz.s = sc.s;
  pil_world_add_overlay(w, sp_overlay, &sc);

  pil_plot *a =
      pil_add_plot(w, 0.62, 0.06, 0.99, 0.46, "t [s]", "\\theta [rad]", 20.0);
  pil_plot_add_signal(a, "\\theta", sig_theta, sc.s);
  pil_plot *b = pil_add_plot(w, 0.62, 0.56, 0.99, 0.97, "t [s]", "", 20.0);
  pil_plot_add_signal(b, "L_z", sig_lz, &lz);
  pil_plot_add_signal(b, "E", sig_E, sc.s);
  pil_plots_record(w, 0.02, 1e-3);

  pil_view view;
  pil_view_fit_region(&view, -1.7, -1.9, 1.7, 0.5, 1300, 760, 0.0, 0.0, 0.6,
                      1.0, 14);
  int rc =
      example_main(w, &view, "Spherical pendulum", argc, argv, "spherical.png");

  pil_world_free(w);
  free(sc.trace);

  return rc;
}
