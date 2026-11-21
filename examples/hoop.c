/*
 * hoop.c - bead on a rotating hoop
 * Hoop of radius R turns about its vertical diameter at constant Omega; bead
 * slides on it. One coordinate, theta from downward vertical
 *   T = 1/2 m R^2 (th'^2 + Omega^2 sin^2 th),   V = -m g R cos th
 * Below Omega^2 = g / R bottom is only equilibrium; above it bottom is unstable
 * and bead oscillates about cos th_0 = g / (R Omega^2), with small-oscillation
 * frequency Omega sin th_0
 * Energy changes (hoop does work) while Jacobi integral h = T_2 - T_0 + V stays
 * constant: both are plotted
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 3.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_lagrange.h>
#include <stdio.h>

typedef struct {
  double m;
  double R;
  double Om;
  double g;
} hoop_params;

typedef struct {
  const hoop_params *p;
  pil_lsys *s;
  const pil_world *w;
  pil_proj3 pr;
} hoop_scene;

static double hoop_T(const double *q, const double *qd, void *user) {
  const hoop_params *p = user;
  double s = sin(q[0]);

  return 0.5 * p->m * p->R * p->R * (qd[0] * qd[0] + p->Om * p->Om * s * s);
}

static double hoop_V(const double *q, void *user) {
  const hoop_params *p = user;

  return -p->m * p->g * p->R * cos(q[0]);
}

// Callback for live angle label
static void theta_label_fn(const pil_world *w, void *user, char *buf, int n) {
  (void)w;
  pil_lsys *s = (pil_lsys *)user;
  double th = pil_lsys_q(s)[0];
  snprintf(buf, n, "\xce\xb8 = %.2f rad", th);
}

// Point of hoop frame (x in hoop plane, z up) turned by phi about z
static pil_vec2 hp(const hoop_scene *sc, double phi, double x, double z) {
  return pil_proj3_point(&sc->pr, x * cos(phi), x * sin(phi), z);
}

static void hoop_overlay(cairo_t *cr, const pil_view *v, void *user) {
  const hoop_scene *sc = user;
  const hoop_params *p = sc->p;
  pil_style st = pil_style_textbook();
  double phi = p->Om * pil_world_time(sc->w);
  double th = pil_lsys_q(sc->s)[0];

  // Axis
  pil_vec2 top =
      pil_view_to_screen(v, pil_proj3_point(&sc->pr, 0, 0, p->R * 1.25));
  pil_vec2 bot =
      pil_view_to_screen(v, pil_proj3_point(&sc->pr, 0, 0, -p->R * 1.25));
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_set_line_width(cr, st.thin_w);
  double dash[2] = {9, 3};
  cairo_set_dash(cr, dash, 2, 0);
  cairo_move_to(cr, top.x, top.y);
  cairo_line_to(cr, bot.x, bot.y);
  cairo_stroke(cr);
  cairo_set_dash(cr, NULL, 0, 0);

  // Hoop
  cairo_set_line_width(cr, st.line_w * 1.4);
  for (int i = 0; i <= 72; i++) {
    double a = 2 * M_PI * i / 72;
    pil_vec2 q =
        pil_view_to_screen(v, hp(sc, phi, p->R * sin(a), -p->R * cos(a)));
    if (i == 0) {
      cairo_move_to(cr, q.x, q.y);
    } else {
      cairo_line_to(cr, q.x, q.y);
    }
  }
  cairo_stroke(cr);

  // Bead and angle against vertical
  pil_vec2 C = pil_view_to_screen(v, pil_proj3_point(&sc->pr, 0, 0, 0));
  pil_vec2 B =
      pil_view_to_screen(v, hp(sc, phi, p->R * sin(th), -p->R * cos(th)));
  
  cairo_set_line_width(cr, st.thin_w);
  cairo_move_to(cr, C.x, C.y);
  cairo_line_to(cr, B.x, B.y);
  cairo_stroke(cr);

  // Angle arc
  double aB = atan2(B.y - C.y, B.x - C.x);
  double aD = atan2(bot.y - C.y, bot.x - C.x);
  double da = aD - aB;
  while (da > M_PI) {
    da -= 2 * M_PI;
  }
  while (da < -M_PI) {
    da += 2 * M_PI;
  }

  cairo_save(cr);
  cairo_set_line_width(cr, st.thin_w);
  cairo_new_path(cr);
  if (da >= 0) {
    cairo_arc(cr, C.x, C.y, 40, aB, aB + da);
  } else {
    cairo_arc(cr, C.x, C.y, 40, aB + da, aB);
  }
  cairo_stroke(cr);
  cairo_restore(cr);

  // Label theta at midpoint of arc
  double midA = aB + da * 0.5;
  pil_draw_text(cr, pil_v2(C.x + 55 * cos(midA), C.y + 55 * sin(midA)),
                "\\theta", 0.5, &st);

  // Draw bead
  cairo_arc(cr, B.x, B.y, 8, 0, 2 * M_PI);
  cairo_set_source_rgb(cr, 1, 1, 1);
  cairo_fill_preserve(cr);
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_set_line_width(cr, st.line_w);
  cairo_stroke(cr);

  pil_draw_text(cr, pil_v2(B.x + 12, B.y - 6), "m", 0.0, &st);
  pil_draw_text(cr, pil_v2(top.x + 12, top.y + 4), "\\Omega", 0.0, &st);
}

static double sig_theta(const pil_world *w, void *u) {
  (void)w;

  return pil_lsys_q(u)[0];
}

static double h0_ref;

static double sig_h(const pil_world *w, void *u) {
  (void)w;

  return pil_lsys_jacobi(u) - h0_ref;
}

typedef struct {
  const hoop_params *p;
  pil_lsys *s;
} energy_ref;

static double E0_ref;

static double sig_E(const pil_world *w, void *u) {
  (void)w;
  const energy_ref *e = u;
  const double *q = pil_lsys_q(e->s), *qd = pil_lsys_qd(e->s);

  return -E0_ref + 0.5 * e->p->m * e->p->R * e->p->R * qd[0] * qd[0] +
         0.5 * e->p->m * e->p->R * e->p->R * e->p->Om * e->p->Om * sin(q[0]) *
             sin(q[0]) -
         e->p->m * e->p->g * e->p->R * cos(q[0]);
}

int main(int argc, char **argv) {
  static hoop_params p = {1.0, 1.0, 5.0, 9.81}; // Omega^2 R = 25 > g
  static hoop_scene sc;
  static energy_ref er;

  pil_world *w = pil_world_new();
  sc.p = &p;
  sc.w = w;
  sc.pr = pil_proj3_make(0, 0, 0.6, 0.35, 1.0);
  sc.s = pil_add_lagrangian(w, 1, hoop_T, hoop_V, &p);
  pil_lsys_set_state(sc.s, (double[]){0.2}, (double[]){0.0});
  er.p = &p;
  er.s = sc.s;
  pil_world_add_overlay(w, hoop_overlay, &sc);

  pil_add_live_label(w, pil_at(-1.2, 1.2), theta_label_fn, sc.s, 0, 0);

  h0_ref = pil_lsys_jacobi(sc.s);
  E0_ref = 0.5 * p.m * p.R * p.R * p.Om * p.Om * sin(0.2) * sin(0.2) -
           p.m * p.g * p.R * cos(0.2);

  pil_plot *a =
      pil_add_plot(w, 0.58, 0.06, 0.99, 0.46, "t [s]", "\\theta [rad]", 6.0);
  pil_plot_add_signal(a, "\\theta", sig_theta, sc.s);
  pil_plot *b =
      pil_add_plot(w, 0.58, 0.56, 0.99, 0.97, "t [s]", "change [J]", 6.0);
  pil_plot_add_signal(b, "E - E_0", sig_E, &er);
  pil_plot_add_signal(b, "h - h_0", sig_h, sc.s);
  pil_plots_record(w, 0.01, 5e-4);

  pil_view view;
  pil_view_fit_region(&view, -1.8, -1.6, 1.8, 1.6, 1300, 760, 0.0, 0.0, 0.55,
                      1.0, 14);
  int rc =
      example_main(w, &view, "Bead on a rotating hoop", argc, argv, "hoop.png");

  pil_world_free(w);

  return rc;
}
