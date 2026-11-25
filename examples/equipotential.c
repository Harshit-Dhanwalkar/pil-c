/*
 * equipotential.c - equipotential surfaces and field lines of an electric
 * dipole Charges +1 and -1 sit on x axis. Thin grey curves are equipotentials
 * V = const, always perpendicular to field lines (black); V = 0 curve is
 * straight line x = 0
 * A probe slides along a horizontal line; thick blue curve is equipotential
 * through it and plot shows V at probe against its position x
 * Moving along an equipotential costs no work; across V changes at rate -E . dl
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 3.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_electrostatic.h>

#define XC 1.2
#define YP 0.9
#define XP 3.6
#define TSW 8.0
#define NMAXP 4000

typedef struct {
  const pil_world *w;
} e_t;

static void state(pil_charges *c) {
  c->n = 2;
  c->q[0] = 1.0;
  c->q[1] = -1.0;
  c->x[0][0] = -XC;
  c->x[0][1] = 0.0;
  c->x[1][0] = XC;
  c->x[1][1] = 0.0;
}

static double probe_x(double t) {
  double s = t / TSW;

  return -XP + 2.0 * XP * (s > 1.0 ? 1.0 : s);
}

static void draw_curve(cairo_t *cr, const pil_view *v, const double *px,
                       const double *py, int n) {
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

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const e_t *e = user;
  pil_style st = pil_style_textbook();
  pil_charges c;
  state(&c);
  static double px[2 * NMAXP + 2], py[2 * NMAXP + 2];
  double xp = probe_x(pil_world_time(e->w));

  // Equipotentials: start points on axis between and beyond charges
  cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
  cairo_set_line_width(cr, st.thin_w);
  double xs[] = {0.0, 0.25, 0.5, 0.8, 1.3, 1.9, 2.7, 3.6, 4.6, 5.4};
  for (unsigned i = 0; i < sizeof xs / sizeof *xs; i++)
    for (int sgn = -1; sgn <= 1; sgn += 2) {
      double x0 = sgn * xs[i];
      if (xs[i] == 0.0 && sgn > 0) {
        continue;
      }
      double y0 = 0.0;
      if (fabs(fabs(x0) - XC) < 0.05) {
        continue;
      }
      double y1 = y0;
      int n =
          pil_es_equipotential(&c, x0, y1, 0.04, NMAXP, 6.0, 3.6, 0.1, px, py);
      draw_curve(cr, v, px, py, n);
    }

  // Equipotential through probe
  cairo_set_source_rgb(cr, 0.1, 0.3, 0.85);
  cairo_set_line_width(cr, 2.4);
  int n = pil_es_equipotential(&c, xp, YP, 0.04, NMAXP, 6.0, 3.6, 0.1, px, py);
  draw_curve(cr, v, px, py, n);

  // Field lines
  cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
  cairo_set_line_width(cr, 1.2);
  for (int j = 0; j < 14; j++) {
    double a = 2 * M_PI * (j + 0.5) / 14;
    n = pil_es_field_line(&c, -XC + 0.12 * cos(a), 0.12 * sin(a), 1, 0, 0.04,
                          NMAXP, 6.0, 3.6, 0.1, px, py);
    draw_curve(cr, v, px, py, n);
  }

  // Charges and probe
  for (int k = 0; k < 2; k++) {
    pil_vec2 q = pil_view_to_screen(v, pil_v2(c.x[k][0], c.x[k][1]));
    if (c.q[k] > 0) {
      cairo_set_source_rgb(cr, 0.85, 0.1, 0.1);
    } else {
      cairo_set_source_rgb(cr, 0.1, 0.25, 0.85);
    }
    cairo_arc(cr, q.x, q.y, 9, 0, 2 * M_PI);
    cairo_fill(cr);
  }

  pil_vec2 pq = pil_view_to_screen(v, pil_v2(xp, YP));
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_arc(cr, pq.x, pq.y, 5, 0, 2 * M_PI);
  cairo_fill(cr);
  cairo_new_path(cr);
}

static double sig_x(const pil_world *w, void *u) {
  (void)u;

  return probe_x(pil_world_time(w));
}

static double sig_v(const pil_world *w, void *u) {
  (void)u;
  pil_charges c;
  state(&c);

  return pil_es_potential(&c, probe_x(pil_world_time(w)), YP);
}

int main(int argc, char **argv) {
  static e_t e;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  e.w = w;
  pil_world_add_overlay(w, overlay, &e);

  pil_plot *p = pil_add_plot(w, 0.69, 0.30, 0.99, 0.70, "x", "V", TSW);
  pil_plot_add_signal(p, "", sig_x, NULL);
  pil_plot_add_signal(p, "V(x, y_p)", sig_v, NULL);
  pil_plot_set_xy(p, 0);
  pil_plots_record(w, 0.02, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -6.0, -3.6, 6.0, 3.6, 1300, 760, 0.0, 0.0, 0.67,
                      1.0, 14);
  int rc = example_main(w, &view, "Equipotentials of a dipole", argc, argv,
                        "equipotential.png");

  pil_world_free(w);

  return rc;
}
