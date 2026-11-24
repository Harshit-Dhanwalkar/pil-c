/*
 * gauss.c - electric field lines of two point charges, and Gauss's law
 * Left charge is +1; right one slowly changes from +1 through 0 to -1
 * Field lines leave positive charges (their number is proportional to charge)
 * and end on negative ones; with equal charges they repel each other, with
 * opposite charges they join up as a dipole field
 * Dashed circle is a Gaussian surface: plot shows flux of E through sphere it
 * represents divided by 4 pi, which equals enclosed charge q_1 + q_2 whatever
 * field lines do inside
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 5.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_electrostatic.h>

#define XC 1.6
#define RG 3.0
#define NLINES 16
#define NMAXP 3000
#define TPER 8.0

typedef struct {
  const pil_world *w;
} g_t;

static void state(double t, pil_charges *c) {
  c->n = 2;
  c->q[0] = 1.0;
  c->q[1] = cos(2.0 * M_PI * t / TPER);
  c->x[0][0] = -XC;
  c->x[0][1] = 0.0;
  c->x[1][0] = XC;
  c->x[1][1] = 0.0;
}

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const g_t *g = user;
  pil_style st = pil_style_textbook();
  pil_charges c;
  state(pil_world_time(g->w), &c);
  static double px[NMAXP + 1], py[NMAXP + 1];

  // Gaussian surface
  double dash[2] = {6, 5};
  pil_vec2 o = pil_view_to_screen(v, pil_v2(0, 0));
  double sc = pil_view_to_screen(v, pil_v2(1, 0)).x - o.x;
  cairo_new_path(cr);
  cairo_set_source_rgb(cr, 0.55, 0.55, 0.55);
  cairo_set_line_width(cr, st.thin_w);
  cairo_set_dash(cr, dash, 2, 0);
  cairo_arc(cr, o.x, o.y, RG * sc, 0, 2 * M_PI);
  cairo_stroke(cr);
  cairo_set_dash(cr, NULL, 0, 0);

  // Field lines: from every charge, number proportional to |q|
  cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
  cairo_set_line_width(cr, 1.3);
  for (int k = 0; k < 2; k++) {
    int nl = (int)lround(NLINES * fabs(c.q[k]));
    int dir = c.q[k] >= 0 ? 1 : -1;
    for (int j = 0; j < nl; j++) {
      double a = 2 * M_PI * (j + 0.5) / nl;
      double x0 = c.x[k][0] + 0.12 * cos(a);
      double y0 = c.x[k][1] + 0.12 * sin(a);
      int n = pil_es_field_line(&c, x0, y0, dir, k, 0.04, NMAXP, 5.5, 3.2, 0.1,
                                px, py);
      for (int i = 0; i < n; i++) {
        pil_vec2 q = pil_view_to_screen(v, pil_v2(px[i], py[i]));
        if (i == 0) {
          cairo_move_to(cr, q.x, q.y);
        } else {
          cairo_line_to(cr, q.x, q.y);
        }
      }

      cairo_stroke(cr);
    }
  }

  // Charges: red positive, blue negative, grey nearly neutral
  for (int k = 0; k < 2; k++) {
    pil_vec2 q = pil_view_to_screen(v, pil_v2(c.x[k][0], c.x[k][1]));
    double a = fmin(1.0, fabs(c.q[k]) * 4);
    if (c.q[k] >= 0) {
      cairo_set_source_rgba(cr, 0.85, 0.1, 0.1, 0.25 + 0.75 * a);
    } else {
      cairo_set_source_rgba(cr, 0.1, 0.25, 0.85, 0.25 + 0.75 * a);
    }
    cairo_arc(cr, q.x, q.y, 6 + 8 * fabs(c.q[k]), 0, 2 * M_PI);
    cairo_fill(cr);
  }

  cairo_new_path(cr);
}

static double sig_flux(const pil_world *w, void *u) {
  (void)u;
  pil_charges c;
  state(pil_world_time(w), &c);

  return pil_es_flux_sphere(&c, 0.0, RG, 200) / (4.0 * M_PI);
}

static double sig_q(const pil_world *w, void *u) {
  (void)u;
  pil_charges c;
  state(pil_world_time(w), &c);

  return c.q[0] + c.q[1];
}

int main(int argc, char **argv) {
  static g_t g;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  g.w = w;
  pil_world_add_overlay(w, overlay, &g);

  pil_plot *p = pil_add_plot(w, 0.66, 0.30, 0.99, 0.70, "t", "charge", TPER);
  pil_plot_add_signal(p, "flux / 4\\pi", sig_flux, &g);
  pil_plot_add_signal(p, "q_1 + q_2", sig_q, &g);
  pil_plots_record(w, 0.05, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -5.5, -3.4, 5.5, 3.4, 1300, 760, 0.0, 0.0, 0.64,
                      1.0, 14);
  int rc = example_main(w, &view, "Field lines and Gauss's law", argc, argv,
                        "gauss.png");

  pil_world_free(w);

  return rc;
}
