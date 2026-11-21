/*
 * doublependulem.c - double pendulum: Lagrangian mechanics and sensitive
 * dependence Two rods, angles theta_1 and theta_2 from vertical Only energies
 * are given to pil_add_lagrangian() T = 1/2 (m_1 + m_2) l_1^2 th1'^2 + 1/2 m_2
 * l_2^2 th2'^2
 *       + m_2 l_1 l_2 th1' th2' cos(th1 - th2)
 *   V = -(m_1 + m_2) g l_1 cos th1 - m_2 g l_2 cos th2
 * A pale twin starts 0.0001 rad away: two agree for a few seconds and then part
 * ways, shown on right as |difference in theta_2|
 * Wake behind lower bob is path of last five seconds
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 14.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_lagrange.h>

#define NTRACE 250
#define TRACE_DT 0.02

typedef struct {
  double m1, m2, l1, l2, g;
} dp_params;

typedef struct {
  const dp_params *p;
  pil_lsys *main_s;
  pil_lsys *twin;
  const pil_world *w;
  double (*trace)[2]; // Lower-bob path, TRACE_DT apart, from t = 0
  int ntrace;
} dp_scene;

static double dp_T(const double *q, const double *qd, void *user) {
  const dp_params *p = user;
  double c = cos(q[0] - q[1]);

  return 0.5 * (p->m1 + p->m2) * p->l1 * p->l1 * qd[0] * qd[0] +
         0.5 * p->m2 * p->l2 * p->l2 * qd[1] * qd[1] +
         p->m2 * p->l1 * p->l2 * qd[0] * qd[1] * c;
}

static double dp_V(const double *q, void *user) {
  const dp_params *p = user;

  return -(p->m1 + p->m2) * p->g * p->l1 * cos(q[0]) -
         p->m2 * p->g * p->l2 * cos(q[1]);
}

static pil_vec2 bob1(const dp_params *p, const double *q) {
  return pil_v2(p->l1 * sin(q[0]), -p->l1 * cos(q[0]));
}

static pil_vec2 bob2(const dp_params *p, const double *q) {
  pil_vec2 a = bob1(p, q);

  return pil_v2(a.x + p->l2 * sin(q[1]), a.y - p->l2 * cos(q[1]));
}

static void draw_pend(cairo_t *cr, const pil_view *v, const dp_params *p,
                      const double *q, double grey, const pil_style *st) {
  pil_vec2 O = pil_view_to_screen(v, pil_v2(0, 0));
  pil_vec2 A = pil_view_to_screen(v, bob1(p, q));
  pil_vec2 B = pil_view_to_screen(v, bob2(p, q));

  cairo_set_source_rgb(cr, grey, grey, grey);
  cairo_set_line_width(cr, st->line_w);
  cairo_move_to(cr, O.x, O.y);
  cairo_line_to(cr, A.x, A.y);
  cairo_line_to(cr, B.x, B.y);
  cairo_stroke(cr);
  pil_vec2 P[2] = {A, B};
  double r[2] = {7 + 3 * p->m1, 7 + 3 * p->m2};
  for (int i = 0; i < 2; i++) {
    cairo_arc(cr, P[i].x, P[i].y, r[i], 0, 2 * M_PI);
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_fill_preserve(cr);
    cairo_set_source_rgb(cr, grey, grey, grey);
    cairo_stroke(cr);
  }
}

static void dp_overlay(cairo_t *cr, const pil_view *v, void *user) {
  const dp_scene *sc = user;
  const dp_params *p = sc->p;
  pil_style st = pil_style_textbook();
  double t = pil_world_time(sc->w);

  // Wake: last NTRACE samples before now, fading with age
  int last = (int)(t / TRACE_DT);
  if (last >= sc->ntrace) {
    last = sc->ntrace - 1;
  }

  int first = last - NTRACE > 0 ? last - NTRACE : 0;
  cairo_set_line_width(cr, 1.1);
  for (int i = first + 1; i <= last; i++) {
    double age = (double)(last - i) / NTRACE, g = 0.15 + 0.8 * age;
    pil_vec2 a =
        pil_view_to_screen(v, pil_v2(sc->trace[i - 1][0], sc->trace[i - 1][1]));
    pil_vec2 b =
        pil_view_to_screen(v, pil_v2(sc->trace[i][0], sc->trace[i][1]));
    cairo_set_source_rgb(cr, g, g, g);
    cairo_move_to(cr, a.x, a.y);
    cairo_line_to(cr, b.x, b.y);
    cairo_stroke(cr);
  }

  draw_pend(cr, v, p, pil_lsys_q(sc->twin), 0.62, &st);
  draw_pend(cr, v, p, pil_lsys_q(sc->main_s), 0.0, &st);

  // Fixed support
  pil_vec2 O = pil_view_to_screen(v, pil_v2(0, 0));
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_set_line_width(cr, st.line_w);
  cairo_move_to(cr, O.x - 22, O.y);
  cairo_line_to(cr, O.x + 22, O.y);
  cairo_stroke(cr);
  for (int i = -2; i <= 2; i++) {
    cairo_move_to(cr, O.x + 9 * i, O.y);
    cairo_line_to(cr, O.x + 9 * i - 6, O.y - 8);
    cairo_stroke(cr);
  }

  pil_draw_text(cr, pil_v2(O.x + 40, O.y + 6), "\\theta_1", 0.0, &st);
}

static double sig_th1(const pil_world *w, void *u) {
  (void)w;

  return pil_lsys_q(u)[0];
}

static double sig_th2(const pil_world *w, void *u) {
  (void)w;

  return pil_lsys_q(u)[1];
}

static double sig_diff(const pil_world *w, void *u) {
  (void)w;
  const dp_scene *sc = u;

  return fabs(pil_lsys_q(sc->main_s)[1] - pil_lsys_q(sc->twin)[1]);
}

int main(int argc, char **argv) {
  static dp_params p = {1.0, 1.0, 1.0, 1.0, 9.81};
  static dp_scene sc;
  const double q0[2] = {2.0, 2.2};

  // Wake precomputed with a world of its own: 5 s of history at any time
  {
    pil_world *w2 = pil_world_new();
    pil_lsys *s2 = pil_add_lagrangian(w2, 2, dp_T, dp_V, &p);
    pil_lsys_set_state(s2, q0, (double[]){0, 0});
    int n = (int)(60.0 / TRACE_DT);
    sc.trace = malloc((size_t)n * sizeof *sc.trace);
    sc.ntrace = n;
    for (int i = 0; i < n; i++) {
      pil_vec2 b = bob2(&p, pil_lsys_q(s2));
      sc.trace[i][0] = b.x;
      sc.trace[i][1] = b.y;
      pil_world_advance(w2, TRACE_DT, 1e-3);
    }

    pil_world_free(w2);
  }

  pil_world *w = pil_world_new();
  sc.p = &p;
  sc.w = w;
  sc.main_s = pil_add_lagrangian(w, 2, dp_T, dp_V, &p);
  sc.twin = pil_add_lagrangian(w, 2, dp_T, dp_V, &p);
  pil_lsys_set_state(sc.main_s, q0, (double[]){0, 0});
  pil_lsys_set_state(sc.twin, (double[]){q0[0] + 1e-4, q0[1]},
                     (double[]){0, 0});
  pil_world_add_overlay(w, dp_overlay, &sc);

  pil_plot *a =
      pil_add_plot(w, 0.58, 0.06, 0.99, 0.46, "t [s]", "\\theta [rad]", 20.0);
  pil_plot_add_signal(a, "\\theta_1", sig_th1, sc.main_s);
  pil_plot_add_signal(a, "\\theta_2", sig_th2, sc.main_s);
  pil_plot *b = pil_add_plot(w, 0.58, 0.56, 0.99, 0.97, "t [s]",
                             "|\\Delta\\theta_2|", 20.0);
  pil_plot_add_signal(b, "twin", sig_diff, &sc);
  pil_plots_record(w, 0.02, 1e-3);

  pil_view view;
  pil_view_fit_region(&view, -2.4, -2.4, 2.4, 0.5, 1300, 800, 0.0, 0.0, 0.55,
                      1.0, 18);
  int rc = example_main(w, &view, "Double pendulum", argc, argv, "dpend.png");

  pil_world_free(w);
  free(sc.trace);

  return rc;
}
