/*
 * rocket.c - rocket equation and staging
 * Two rockets of same total mass (1000) and same payload (50), exhaust speed ve
 * = 2500 and mdot = 10 per second, rise in uniform gravity g = 9.81
 * Left: one stage, 150 of structure and 800 of propellant
 * Right: two stages (100 + 600 and 50 + 200); empty first stage is dropped at t
 * = 60
 * Dropping dead weight pays: without gravity velocity gain is ve ln(m_0 / m_1)
 * per stage and two-stage rocket ends faster and higher
 * Plots: velocity and altitude of both
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 20.0
#define PIL_EXAMPLE_IMAGE_TIME 70.0
#include <math.h>
#include <pil/pil_rocket.h>

#include "example_main.h"

#define TEND 300.0
#define DT 0.1

static const pil_rocket R[2] = {
    {9.81, 50.0, 1, {{150.0, 800.0, 10.0, 2500.0}}},
    {9.81,
     50.0,
     2,
     {{100.0, 600.0, 10.0, 2500.0}, {50.0, 200.0, 10.0, 2500.0}}}};
static pil_rocket_state S[2];
static double tnow = -1.0;

static void advance(const pil_world *w) {
  double t = pil_world_time(w);
  if (tnow < 0.0 || t < tnow) {
    pil_rocket_init(&S[0]);
    pil_rocket_init(&S[1]);
    tnow = 0.0;
  }

  while (t - tnow > 1e-9) {
    double h = t - tnow < DT ? t - tnow : DT;
    pil_rocket_step(&R[0], &S[0], h);
    pil_rocket_step(&R[1], &S[1], h);
    tnow += h;
  }
}

typedef struct {
  const pil_world *w;
} r_t;

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const r_t *g = user;
  pil_style st = pil_style_textbook();
  advance(g->w);
  cairo_new_path(cr);
  cairo_set_line_width(cr, st.thin_w);
  cairo_set_source_rgb(cr, 0.7, 0.7, 0.7);
  for (int k = 0; k <= 8; k++) {
    pil_vec2 a = pil_view_to_screen(v, pil_v2(-3.6, k));
    pil_vec2 b = pil_view_to_screen(v, pil_v2(3.6, k));
    cairo_move_to(cr, a.x, a.y);
    cairo_line_to(cr, b.x, b.y);
    cairo_stroke(cr);
    char lab[24];
    snprintf(lab, sizeof lab, "%d km", k * 100);
    cairo_set_source_rgb(cr, 0.35, 0.35, 0.35);
    pil_draw_text(cr, pil_view_to_screen(v, pil_v2(-3.7, k + 0.12)), lab, 0.0,
                  &st);
    cairo_set_source_rgb(cr, 0.7, 0.7, 0.7);
  }

  for (int i = 0; i < 2; i++) {
    double x = i == 0 ? -1.6 : 1.6;
    double h = S[i].h / 1e5;
    int burning = pil_rocket_burning(&R[i], &S[i]);
    // Body (taller while more stages are attached)
    double body =
        0.16 + 0.06 * (R[i].n - S[i].stage > 0 ? R[i].n - S[i].stage : 0);
    pil_vec2 p0 = pil_view_to_screen(v, pil_v2(x - 0.12, h));
    pil_vec2 p1 = pil_view_to_screen(v, pil_v2(x + 0.12, h));
    pil_vec2 p2 = pil_view_to_screen(v, pil_v2(x, h + body));
    cairo_new_path(cr);
    cairo_move_to(cr, p0.x, p0.y);
    cairo_line_to(cr, p1.x, p1.y);
    cairo_line_to(cr, p2.x, p2.y);
    cairo_close_path(cr);
    cairo_set_source_rgb(cr, 0.2, 0.3, 0.7);
    cairo_fill(cr);
    if (burning) {
      pil_vec2 f = pil_view_to_screen(v, pil_v2(x, h - 0.28));
      cairo_new_path(cr);
      cairo_move_to(cr, p0.x, p0.y);
      cairo_line_to(cr, p1.x, p1.y);
      cairo_line_to(cr, f.x, f.y);
      cairo_close_path(cr);
      cairo_set_source_rgb(cr, 1.0, 0.55, 0.1);
      cairo_fill(cr);
    }

    char buf[48];
    snprintf(buf, sizeof buf, i == 0 ? "1 stage" : "2 stages (stage %d)",
             S[i].stage < R[i].n ? S[i].stage + 1 : R[i].n);
    cairo_set_source_rgb(cr, 0, 0, 0);
    pil_draw_text(cr, pil_view_to_screen(v, pil_v2(x, -0.55)), buf, 0.5, &st);
  }

  cairo_new_path(cr);
}

static double s_v(const pil_world *w, void *u) {
  advance(w);

  return S[*(int *)u].v / 1e3;
}

static double s_h(const pil_world *w, void *u) {
  advance(w);

  return S[*(int *)u].h / 1e3;
}

int main(int argc, char **argv) {
  static r_t g;
  static int i0 = 0;
  static int i1 = 1;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  g.w = w;
  pil_world_add_overlay(w, overlay, &g);

  pil_plot *p =
      pil_add_plot(w, 0.62, 0.06, 0.99, 0.47, "t [s]", "v [km/s]", TEND);
  pil_plot_add_signal(p, "1 stage", s_v, &i0);
  pil_plot_add_signal(p, "2 stages", s_v, &i1);
  pil_plot *q =
      pil_add_plot(w, 0.62, 0.56, 0.99, 0.97, "t [s]", "h [km]", TEND);
  pil_plot_add_signal(q, "1 stage", s_h, &i0);
  pil_plot_add_signal(q, "2 stages", s_h, &i1);
  pil_plots_record(w, 1.0, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -4.4, -1.0, 3.6, 8.6, 1300, 760, 0.0, 0.0, 0.6,
                      1.0, 14);
  int rc = example_main(w, &view, "Rocket equation and staging", argc, argv,
                        "rocket.png");

  pil_world_free(w);

  return rc;
}
