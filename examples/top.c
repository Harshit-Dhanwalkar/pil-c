/*
 * top.c - heavy symmetric top on a fixed point, in Euler angles
 * A flywheel spins fast (psi) on an axis that leans by theta from vertical and
 * turns about it (phi); a hand-held look at precession and nutation. Problem is
 * 3D, so it is stated as kinetic and potential energy and handed to
 * pil_add_lagrangian(); two cyclic angles give conserved momenta p_psi and
 * p_phi, written under figure with one remaining equation
 * Everything is public-API
 */

#include "example_main.h"
#include <math.h>
#include <pil/pil_lagrange.h>

typedef struct {
  double I1, I3, m, l, g;
} top_params;

typedef struct {
  const top_params *p;
  pil_lsys *s;
  pil_proj3 pr;
  double axis, radius;
} top_scene;

// q = (phi, theta, psi): turn about vertical, lean from vertical, spin
static double top_T(const double *q, const double *qd, void *user) {
  const top_params *p = user;
  double st = sin(q[1]);
  double w3 = qd[2] + qd[0] * cos(q[1]);

  return 0.5 * p->I1 * (qd[1] * qd[1] + qd[0] * qd[0] * st * st) +
         0.5 * p->I3 * w3 * w3;
}

static double top_V(const double *q, void *user) {
  const top_params *p = user;

  return p->m * p->g * p->l * cos(q[1]);
}

static double sig_theta(const pil_world *w, void *user) {
  (void)w;

  return pil_lsys_q(user)[1];
}

static double sig_phidot(const pil_world *w, void *user) {
  (void)w;

  return pil_lsys_qd(user)[0];
}

static pil_vec2 scr(const top_scene *sc, const pil_view *v, const double *a,
                    double s) {
  return pil_view_to_screen(
      v, pil_proj3_point(&sc->pr, s * a[0], s * a[1], s * a[2]));
}

static void arc_px(cairo_t *cr, const top_scene *sc, const pil_view *v,
                   double r, double z, double a0, double a1) {
  for (int i = 0; i <= 32; i++) {
    double a = a0 + (a1 - a0) * i / 32.0;
    double p[3] = {r * cos(a), r * sin(a), z};
    pil_vec2 s = scr(sc, v, p, 1.0);
    if (i == 0) {
      cairo_move_to(cr, s.x, s.y);
    } else {
      cairo_line_to(cr, s.x, s.y);
    }
  }
  cairo_stroke(cr);
}

static void top_overlay(cairo_t *cr, const pil_view *v, void *user) {
  const top_scene *sc = user;
  const top_params *p = sc->p;
  pil_style st = pil_style_textbook();
  const double *q = pil_lsys_q(sc->s);
  double ph = q[0];
  double th = q[1];
  double ps = q[2];
  // Figure axis and two directions perpendicular to it; e1 is line of nodes, e2
  // completes frame
  double a[3] = {sin(th) * sin(ph), -sin(th) * cos(ph), cos(th)};
  double e1[3] = {cos(ph), sin(ph), 0};
  double e2[3] = {a[1] * e1[2] - a[2] * e1[1], a[2] * e1[0] - a[0] * e1[2],
                  a[0] * e1[1] - a[1] * e1[0]};
  double o[3] = {0, 0, 0};
  double L = sc->axis;
  double R = sc->radius;
  double cen[3] = {a[0] * 0.8 * L, a[1] * 0.8 * L, a[2] * 0.8 * L};
  double tip[3] = {a[0] * L, a[1] * L, a[2] * L};
  pil_vec2 O = scr(sc, v, o, 1.0);
  pil_vec2 T = scr(sc, v, tip, 1.0);
  double dash[] = {6, 4};

  cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);

  // Floor, vertical axis, and angle phi swept in horizontal plane
  cairo_set_line_width(cr, st.thin_w);
  arc_px(cr, sc, v, 1.25 * R, -0.0, 0, 2 * M_PI);
  cairo_set_dash(cr, dash, 2, 0);
  {
    double top[3] = {0, 0, 1.15 * L};
    pil_vec2 Z = scr(sc, v, top, 1.0);
    cairo_move_to(cr, O.x, O.y);
    cairo_line_to(cr, Z.x, Z.y);
    cairo_stroke(cr);
    // line of nodes
    double n[3] = {1.25 * R * e1[0], 1.25 * R * e1[1], 0};
    pil_vec2 N = scr(sc, v, n, 1.0);
    cairo_move_to(cr, O.x, O.y);
    cairo_line_to(cr, N.x, N.y);
    cairo_stroke(cr);
    // lean plane: drop from tip to horizontal plane
    double foot[3] = {tip[0], tip[1], 0};
    pil_vec2 Fp = scr(sc, v, foot, 1.0);
    cairo_move_to(cr, T.x, T.y);
    cairo_line_to(cr, Fp.x, Fp.y);
    cairo_line_to(cr, O.x, O.y);
    cairo_stroke(cr);
    pil_draw_text(cr, (pil_vec2){Z.x + 8, Z.y + 4}, "z", 0, &st);
  }
  cairo_set_dash(cr, NULL, 0, 0);

  // Angle markers: phi about vertical (from x axis), theta from vertical
  arc_px(cr, sc, v, 0.45 * R, 0, 0, ph + M_PI / 2 > 0 ? ph : ph);
  {
    double m[3] = {0.55 * R * cos(ph / 2), 0.55 * R * sin(ph / 2), 0};
    pil_vec2 s = scr(sc, v, m, 1.0);
    pil_draw_text(cr, (pil_vec2){s.x + 6, s.y + 14}, "\xcf\x86", 0, &st);
  }
  {
    double m[3] = {0.5 * L * sin(th / 2) * sin(ph),
                   -0.5 * L * sin(th / 2) * cos(ph), 0.5 * L * cos(th / 2)};
    pil_vec2 s = scr(sc, v, m, 1.0);
    pil_draw_text(cr, (pil_vec2){s.x - 8, s.y - 6}, "\xce\xb8", 1, &st);
  }

  // Rod, then flywheel as a polygon in plane (e1, e2) around cen; rod goes on
  // top when wheel is behind it
  int rod_front = pil_proj3_depth(&sc->pr, cen[0], cen[1], cen[2]) <
                          pil_proj3_depth(&sc->pr, 0, 0, 0) + 1e-9
                      ? 0
                      : 1;
  cairo_set_line_width(cr, 3.0);
  if (!rod_front) {
    cairo_move_to(cr, O.x, O.y);
    cairo_line_to(cr, T.x, T.y);
    cairo_stroke(cr);
  }
  cairo_set_line_width(cr, st.line_w);
  cairo_set_line_width(cr, st.line_w);
  for (int i = 0; i < 64; i++) {
    double t = 2 * M_PI * i / 64.0;
    double c = cos(t);
    double s = sin(t);
    double pt[3];
    for (int k = 0; k < 3; k++) {
      pt[k] = cen[k] + R * (c * e1[k] + s * e2[k]);
    }
    pil_vec2 sp = scr(sc, v, pt, 1.0);
    if (i == 0) {
      cairo_move_to(cr, sp.x, sp.y);
    } else {
      cairo_line_to(cr, sp.x, sp.y);
    }
  }
  cairo_close_path(cr);
  cairo_set_source_rgb(cr, st.body_fill, st.body_fill, st.body_fill);
  cairo_fill_preserve(cr);
  cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);
  cairo_stroke(cr);

  // Spokes turning with spin angle psi
  cairo_set_line_width(cr, st.thin_w);
  for (int i = 0; i < 4; i++) {
    double t = ps + i * M_PI / 2;
    double pt[3];
    for (int k = 0; k < 3; k++) {
      pt[k] = cen[k] + R * (cos(t) * e1[k] + sin(t) * e2[k]);
    }
    pil_vec2 c0 = scr(sc, v, cen, 1.0);
    pil_vec2 sp = scr(sc, v, pt, 1.0);
    cairo_move_to(cr, c0.x, c0.y);
    cairo_line_to(cr, sp.x, sp.y);
    cairo_stroke(cr);
  }
  if (rod_front) {
    cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);
    cairo_set_line_width(cr, 3.0);
    cairo_move_to(cr, O.x, O.y);
    cairo_line_to(cr, T.x, T.y);
    cairo_stroke(cr);
  }
  // pivot
  pil_draw_hinge(cr, O, &st);
  {
    double lb[3] = {cen[0] + 1.15 * R * e1[0], cen[1] + 1.15 * R * e1[1],
                    cen[2] + 1.15 * R * e1[2]};
    pil_vec2 s = scr(sc, v, lb, 1.0);
    pil_draw_text(cr, (pil_vec2){s.x + 6, s.y}, "m", 0, &st);
    pil_vec2 c0 = scr(sc, v, cen, 1.0);
    (void)c0;
  }
  {
    double lb[3] = {a[0] * 0.3 * L, a[1] * 0.3 * L, a[2] * 0.3 * L};
    pil_vec2 s = scr(sc, v, lb, 1.0);
    pil_draw_text(cr, (pil_vec2){s.x + 10, s.y + 4}, "l", 0, &st);
  }
  (void)p;

  // Equations under figure
  double H = v->height;
  double y0 = H * 0.70;
  cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);
  pil_draw_text(cr, (pil_vec2){24, y0},
                "L = 1/2 I_1 (\xce\xb8'^2 + \xcf\x86'^2 sin^2 \xce\xb8) + "
                "1/2 I_3 (\xcf\x88' + \xcf\x86' cos \xce\xb8)^2 - m g l cos "
                "\xce\xb8",
                0, &st);
  pil_draw_text(cr, (pil_vec2){24, y0 + 40},
                "p_\xcf\x88 = I_3 (\xcf\x88' + \xcf\x86' cos \xce\xb8) = const"
                "     (\xcf\x88 is cyclic)",
                0, &st);
  pil_draw_text(cr, (pil_vec2){24, y0 + 80},
                "p_\xcf\x86 = I_1 \xcf\x86' sin^2 \xce\xb8 + p_\xcf\x88 cos "
                "\xce\xb8 = const     (\xcf\x86 is cyclic)",
                0, &st);
  pil_draw_text(cr, (pil_vec2){24, y0 + 120},
                "I_1 \xce\xb8'' = I_1 \xcf\x86'^2 sin \xce\xb8 cos \xce\xb8 - "
                "p_\xcf\x88 \xcf\x86' sin \xce\xb8 + m g l sin \xce\xb8",
                0, &st);
  pil_draw_text(cr, (pil_vec2){24, y0 + 165},
                "fast spin: \xcf\x86' ~ m g l / p_\xcf\x88"
                " (precession), small wobble in \xce\xb8 (nutation)",
                0, &st);
}

int main(int argc, char **argv) {
  static top_params p = {0.03, 0.01, 1.0, 0.3, 9.81};
  static top_scene sc;
  pil_world *w = pil_world_new();

  pil_lsys *s = pil_add_lagrangian(w, 3, top_T, top_V, &p);
  pil_lsys_set_state(s, (double[]){2.3, 0.7, 0.0}, (double[]){4.2, 0.0, 70.0});

  sc.p = &p;
  sc.s = s;
  sc.axis = 1.7;
  sc.radius = 0.5;
  sc.pr = pil_proj3_make(0, 0, 0.5, 0.45, 1.0);
  pil_world_add_overlay(w, top_overlay, &sc);

  pil_plot *pl = pil_add_plot(w, 0.66, 0.04, 0.99, 0.42, "t [s]",
                              "\xce\xb8 [rad], \xcf\x86' [rad/s]", 4.0);
  pil_plot_add_signal(pl, "\xce\xb8", sig_theta, s);
  pil_plot_add_signal(pl, "\xcf\x86'", sig_phidot, s);
  pil_plots_record(w, 0.01, 5e-4);

  pil_view view;
  pil_view_fit_region(&view, -1.7, -0.8, 1.7, 1.9, 1200, 1000, 0.0, 0.0, 0.62,
                      0.66, 6);
  int rc = example_main(w, &view, "Heavy symmetric top", argc, argv, "top.png");

  pil_world_free(w);

  return rc;
}
