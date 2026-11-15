/*
 * tower.c - rotating tower: a block on a radial arm and a pendulum under it
 * Arm turns about vertical axis at constant Omega
 * Block m_1 (coordinate x) slides on arm against spring k; pendulum m_2 (length
 * l, angle gamma) hangs from block and swings across arm
 * Problem is 3D, so it s stated as energies T and V and handed to
 * pil_add_lagrangian(); Coriolis and centrifugal terms come out of
 * Euler-Lagrange equations, shown under figure as M q'' + G q' + F(q) = 0
 * Everything is public-API
 */

#include "example_main.h"
#include <math.h>
#include <pil/pil_lagrange.h>
#include <string.h>

typedef struct {
  double m1, m2, l, k, L0, Om, g;
} tower_params;

typedef struct {
  const tower_params *p;
  pil_lsys *s;
  pil_proj3 pr;
  const pil_world *w;
  double height, arm;
} tower_scene;

// Block distance x along arm and pendulum angle gamma; arm turns at Om
static double tower_T(const double *q, const double *qd, void *user) {
  const tower_params *p = user;
  double x = q[0];
  double xd = qd[0];
  double gd = qd[1];
  double sg = sin(q[1]);
  double cg = cos(q[1]);
  double vx = xd - p->Om * p->l * sg;
  double vy = p->l * cg * gd + p->Om * x;
  double vz = p->l * sg * gd;

  return 0.5 * p->m1 * (xd * xd + p->Om * p->Om * x * x) +
         0.5 * p->m2 * (vx * vx + vy * vy + vz * vz);
}

static double tower_V(const double *q, void *user) {
  const tower_params *p = user;

  return 0.5 * p->k * (q[0] - p->L0) * (q[0] - p->L0) -
         p->m2 * p->g * p->l * cos(q[1]);
}

static double sig_x(const pil_world *w, void *user) {
  (void)w;

  return pil_lsys_q(user)[0];
}

static double sig_gamma(const pil_world *w, void *user) {
  (void)w;

  return pil_lsys_q(user)[1];
}

// Point of rotating frame (x along arm, y across, z up) seen from outside
static pil_vec2 scene_point(const tower_scene *sc, double th, double x,
                            double y, double z) {
  double c = cos(th);
  double s = sin(th);

  return pil_proj3_point(&sc->pr, x * c - y * s, x * s + y * c, z);
}

static void line_px(cairo_t *cr, pil_vec2 a, pil_vec2 b) {
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, b.x, b.y);
  cairo_stroke(cr);
}

static void ellipse_ring(cairo_t *cr, const tower_scene *sc, const pil_view *v,
                         double r, double z, double a0, double a1) {
  for (int i = 0; i <= 48; i++) {
    double a = a0 + (a1 - a0) * i / 48.0;
    pil_vec2 p = pil_view_to_screen(
        v, pil_proj3_point(&sc->pr, r * cos(a), r * sin(a), z));
    if (i == 0) {
      cairo_move_to(cr, p.x, p.y);
    } else {
      cairo_line_to(cr, p.x, p.y);
    }
  }
  cairo_stroke(cr);
}

// A matrix of text cells, brackets on both sides; returns right edge
static double matrix_px(cairo_t *cr, double x, double yc, int rows, int cols,
                        const char *const *cell, const double *colw,
                        const pil_style *st) {
  double rowh = st->font_px * 1.9;
  double w = 0;
  double top = yc - rows * rowh / 2;
  double bot = yc + rows * rowh / 2;
  for (int c = 0; c < cols; c++) {
    w += colw[c];
  }

  double x0 = x + 8;
  double x1 = x0 + w + 8;
  cairo_set_line_width(cr, st->thin_w);
  cairo_move_to(cr, x0 - 4, top);
  cairo_line_to(cr, x, top);
  cairo_line_to(cr, x, bot);
  cairo_line_to(cr, x0 - 4, bot);
  cairo_move_to(cr, x1 + 4, top);
  cairo_line_to(cr, x1 + 12, top);
  cairo_line_to(cr, x1 + 12, bot);
  cairo_line_to(cr, x1 + 4, bot);
  cairo_stroke(cr);

  double cx = x0 + 4;
  for (int c = 0; c < cols; c++) {
    for (int r = 0; r < rows; r++) {
      pil_draw_text(
          cr,
          (pil_vec2){cx + colw[c] / 2, yc + (r - (rows - 1) / 2.0) * rowh + 6},
          cell[r * cols + c], 0.5, st);
    }
    cx += colw[c];
  }

  return x1 + 12;
}

static void tower_overlay(cairo_t *cr, const pil_view *v, void *user) {
  const tower_scene *sc = user;
  const tower_params *p = sc->p;
  pil_style st = pil_style_textbook();
  const double *q = pil_lsys_q(sc->s);
  double th = p->Om * pil_world_time(sc->w);
  double x = q[0];
  double gm = q[1];
  double yb = p->l * sin(gm);
  double zb = -p->l * cos(gm);
  pil_vec2 O = pil_view_to_screen(v, pil_proj3_point(&sc->pr, 0, 0, 0));
  pil_vec2 base =
      pil_view_to_screen(v, pil_proj3_point(&sc->pr, 0, 0, -sc->height));
  pil_vec2 blk = pil_view_to_screen(v, scene_point(sc, th, x, 0, 0));
  pil_vec2 tip = pil_view_to_screen(v, scene_point(sc, th, sc->arm, 0, 0));
  pil_vec2 bob = pil_view_to_screen(v, scene_point(sc, th, x, yb, zb));

  cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);
  cairo_set_line_width(cr, st.line_w);

  // Ground ring and tower
  cairo_set_line_width(cr, st.thin_w);
  ellipse_ring(cr, sc, v, 0.55 * sc->arm, -sc->height, 0, 2 * M_PI);
  cairo_set_line_width(cr, st.line_w);
  cairo_set_source_rgb(cr, 0.82, 0.82, 0.82);
  cairo_rectangle(cr, O.x - 7, O.y, 14, base.y - O.y);
  cairo_fill_preserve(cr);
  cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);
  cairo_stroke(cr);

  // Arm
  cairo_set_line_width(cr, 5.0);
  line_px(cr, O, tip);
  cairo_set_line_width(cr, st.line_w);

  // Axis, rotation arrow Omega about it
  cairo_set_line_width(cr, st.thin_w);
  double dash[] = {6, 4};
  cairo_set_dash(cr, dash, 2, 0);
  line_px(cr, (pil_vec2){O.x, O.y - 40}, O);
  cairo_set_dash(cr, NULL, 0, 0);
  ellipse_ring(cr, sc, v, 0.4 * sc->arm, 0.5, 0.5, 5.0);
  {
    double a = 5.0;
    double r = 0.4 * sc->arm;
    pil_vec2 e = pil_view_to_screen(
        v, pil_proj3_point(&sc->pr, r * cos(a), r * sin(a), 0.5));
    pil_vec2 b = pil_view_to_screen(
        v, pil_proj3_point(&sc->pr, r * cos(a - 0.12), r * sin(a - 0.12), 0.5));
    pil_draw_arrow(cr, b, e, 9, &st);
    pil_draw_text(cr, (pil_vec2){e.x + 12, e.y - 4}, "\xce\xa9", 0, &st);
  }

  // Spring from axis to block
  cairo_set_line_width(cr, st.thin_w);
  pil_draw_spring(cr, O, blk, 9, 7, &st);
  {
    pil_vec2 mid = {(O.x + blk.x) / 2, (O.y + blk.y) / 2};
    pil_draw_text(cr, (pil_vec2){mid.x - 4, mid.y + 24}, "k", 0.5, &st);
    pil_draw_text(cr, (pil_vec2){(blk.x + tip.x) / 2, (blk.y + tip.y) / 2 + 22},
                  "x", 0.5, &st);
  }

  // Block, rod, bob
  cairo_set_line_width(cr, st.line_w);
  cairo_set_source_rgb(cr, st.body_fill, st.body_fill, st.body_fill);
  cairo_rectangle(cr, blk.x - 11, blk.y - 11, 22, 22);
  cairo_fill_preserve(cr);
  cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);
  cairo_stroke(cr);
  pil_draw_text(cr, (pil_vec2){blk.x - 8, blk.y - 20}, "m_1", 0, &st);

  pil_vec2 bobb = bob;
  cairo_set_line_width(cr, 2.0);
  line_px(cr, blk, bobb);
  cairo_set_line_width(cr, st.line_w);
  cairo_set_source_rgb(cr, st.body_fill, st.body_fill, st.body_fill);
  cairo_arc(cr, bobb.x, bobb.y, 10, 0, 2 * M_PI);
  cairo_fill_preserve(cr);
  cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);
  cairo_stroke(cr);
  pil_draw_text(cr, (pil_vec2){bobb.x + 15, bobb.y + 6}, "m_2", 0, &st);
  pil_draw_text(cr, (pil_vec2){(blk.x + bobb.x) / 2 - 12, (blk.y + bobb.y) / 2},
                "l", 1, &st);

  // Angle gamma between vertical through block and rod
  {
    pil_vec2 down = pil_view_to_screen(v, scene_point(sc, th, x, 0, -p->l));
    cairo_set_line_width(cr, st.thin_w);
    cairo_set_dash(cr, dash, 2, 0);
    line_px(cr, blk, down);
    cairo_set_dash(cr, NULL, 0, 0);
    pil_draw_text(cr, (pil_vec2){down.x + 8, down.y + 18}, "\xce\xb3", 0, &st);
  }

  // Equations of motion under figure
  double W = v->width;
  double H = v->height;
  double y0 = H * 0.68;
  cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);
  pil_draw_text(cr, (pil_vec2){24, y0 - 34},
                "M q'' + G q' + F(q) = 0,    q = (x, \xce\xb3)", 0, &st);
  const char *M[] = {"m_1 + m_2", "0", "0", "m_2 l^2"};
  const char *G[] = {"0", "-2 m_2 \xce\xa9 l cos \xce\xb3",
                     "2 m_2 \xce\xa9 l cos \xce\xb3", "0"};
  const char *F[] = {"k (x - L_0) - (m_1 + m_2) \xce\xa9^2 x",
                     "m_2 l (g - \xce\xa9^2 l cos \xce\xb3) sin \xce\xb3"};
  const char *qdd[] = {"x''", "\xce\xb3''"};
  const char *qd[] = {"x'", "\xce\xb3'"};
  const char *zero[] = {"0", "0"};
  double cwM[] = {84, 64};
  double cwG[] = {150, 150};
  double cwq[] = {40};
  double cwF[] = {290};
  double cwz[] = {14};
  double xe = 24;
  (void)W;
  xe = matrix_px(cr, xe, y0 + 20, 2, 2, M, cwM, &st);
  xe = matrix_px(cr, xe + 6, y0 + 20, 2, 1, qdd, cwq, &st);
  pil_draw_text(cr, (pil_vec2){xe + 22, y0 + 26}, "+", 0.5, &st);
  xe = matrix_px(cr, xe + 44, y0 + 20, 2, 2, G, cwG, &st);
  xe = matrix_px(cr, xe + 6, y0 + 20, 2, 1, qd, cwq, &st);
  pil_draw_text(cr, (pil_vec2){24, y0 + 100}, "+", 0, &st);
  xe = matrix_px(cr, 60, y0 + 120, 2, 1, F, cwF, &st);
  pil_draw_text(cr, (pil_vec2){xe + 20, y0 + 126}, "=", 0.5, &st);
  matrix_px(cr, xe + 40, y0 + 120, 2, 1, zero, cwz, &st);
  pil_draw_text(cr, (pil_vec2){24, y0 + 190},
                "G is skew: Coriolis forces do no work; "
                "-(m_1 + m_2) \xce\xa9^2 x is centrifugal push",
                0, &st);
}

int main(int argc, char **argv) {
  static tower_params p = {1.0, 0.6, 0.5, 60.0, 0.8, 2.0, 9.81};
  static tower_scene sc;
  pil_world *w = pil_world_new();

  pil_lsys *s = pil_add_lagrangian(w, 2, tower_T, tower_V, &p);
  pil_lsys_set_state(s, (double[]){0.9, 0.7}, (double[]){0.3, -0.5});

  sc.p = &p;
  sc.s = s;
  sc.w = w;
  sc.height = 1.1;
  sc.arm = 1.5;
  sc.pr = pil_proj3_make(0, 0, 0.5, 0.45, 1.0);
  pil_world_add_overlay(w, tower_overlay, &sc);

  pil_plot *pl = pil_add_plot(w, 0.66, 0.04, 0.99, 0.42, "t [s]",
                              "x [m], \xce\xb3 [rad]", 8.0);
  pil_plot_add_signal(pl, "x", sig_x, s);
  pil_plot_add_signal(pl, "\xce\xb3", sig_gamma, s);
  pil_plots_record(w, 0.05, 1e-3);

  pil_view view;
  pil_view_fit_region(&view, -1.7, -1.5, 1.7, 1.0, 1200, 1000, 0.0, 0.0, 0.62,
                      0.66, 6);
  int rc = example_main(w, &view, "Rotating tower with a sliding pendulum",
                        argc, argv, "tower.png");

  pil_world_free(w);

  return rc;
}
