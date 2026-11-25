/*
 * prism.c - dispersion of white light by a glass prism
 * A narrow beam of white light hits an equilateral prism (apex 60 degrees,
 * dense flint glass: Cauchy n = a + b / lambda^2). Refraction at two
 * faces bends short wavelengths more than long ones, so beam fans out
 * into a spectrum, violet deviated most. Angle of incidence is swept;
 * plot shows deviation of violet (400 nm) and red (700 nm) light, each
 * with its minimum where ray crosses prism symmetrically
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 3.4
#include "example_main.h"
#include <math.h>
#include <pil/pil_optics.h>

#define APEX (M_PI / 3.0)
#define CA 1.60
#define CB 0.010
#define TSWEEP 10.0
#define NRAY 13

typedef struct {
  const pil_world *w;
} pr_t;

static double inc(double t) {
  double s = t / TSWEEP;

  return (45.0 + 30.0 * (s > 1 ? 1 : s)) * M_PI / 180.0; // 45 to 75 degrees
}

static double idx(double nm) { return pil_cauchy_index(CA, CB, nm * 1e-3); }

static const double TRI[3][2] = {
    {-1.5, 0.0}, {1.5, 0.0}, {0.0, 2.598076}}; // CCW

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const pr_t *p = user;
  pil_style st = pil_style_textbook();
  double ti = inc(pil_world_time(p->w));

  // prism
  cairo_new_path(cr);
  for (int i = 0; i < 3; i++) {
    pil_vec2 q = pil_view_to_screen(v, pil_v2(TRI[i][0], TRI[i][1]));
    if (i == 0)
      cairo_move_to(cr, q.x, q.y);
    else
      cairo_line_to(cr, q.x, q.y);
  }

  cairo_close_path(cr);
  cairo_set_source_rgba(cr, 0.6, 0.75, 0.95, 0.3);
  cairo_fill_preserve(cr);
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_set_line_width(cr, st.line_w);
  cairo_stroke(cr);

  // Entry point on left face (apex -> base-left edge: from vertex 2 to vertex
  // 0), middle
  double pe[2] = {0.5 * (TRI[0][0] + TRI[2][0]), 0.5 * (TRI[0][1] + TRI[2][1])};
  // Incident direction: inward normal of left face rotated by incidence angle
  double inward = -30.0 * M_PI / 180.0, ang = inward + ti;
  double d[2] = {cos(ang), sin(ang)},
         start[2] = {pe[0] - 4.2 * d[0], pe[1] - 4.2 * d[1]};

  // White beam, then each colour
  static double pts[8][2];
  cairo_set_line_width(cr, 3.0);
  cairo_set_source_rgb(cr, 0.15, 0.15, 0.15);
  pil_vec2 a = pil_view_to_screen(v, pil_v2(start[0], start[1]));
  pil_vec2 b = pil_view_to_screen(v, pil_v2(pe[0], pe[1]));
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, b.x, b.y);
  cairo_stroke(cr);
  for (int k = 0; k < NRAY; k++) {
    double nm = 400.0 + 300.0 * k / (NRAY - 1.0);
    double rgb[3];
    pil_wavelength_rgb(nm, rgb);
    int np = pil_glass_trace(TRI, 3, idx(nm), start, d, 6.0, pts, 8, NULL);
    cairo_set_source_rgb(cr, rgb[0] * 0.9, rgb[1] * 0.9, rgb[2] * 0.9);
    cairo_set_line_width(cr, 1.8);
    for (int i = 1; i < np; i++) { // Skip part before glass (shared white beam)
      pil_vec2 q = pil_view_to_screen(v, pil_v2(pts[i][0], pts[i][1]));
      if (i == 1) {
        pil_vec2 q0 = pil_view_to_screen(v, pil_v2(pts[0][0], pts[0][1]));
        (void)q0;
        cairo_move_to(cr, q.x, q.y);
      } else {
        cairo_line_to(cr, q.x, q.y);
      }
    }

    cairo_stroke(cr);
  }

  char buf[64];
  snprintf(buf, sizeof buf, "\\theta_i = %.1f", ti * 180.0 / M_PI);
  pil_draw_text(cr, pil_view_to_screen(v, pil_v2(-6.5, 2.7)), buf, 0.0, &st);
  cairo_new_path(cr);
}

static double sig_deg(const pil_world *w, void *u) {
  (void)u;

  return inc(pil_world_time(w)) * 180.0 / M_PI;
}

static double sig_v(const pil_world *w, void *u) {
  (void)u;

  return pil_prism_deviation(APEX, idx(400.0), inc(pil_world_time(w))) * 180.0 /
         M_PI;
}

static double sig_r(const pil_world *w, void *u) {
  (void)u;

  return pil_prism_deviation(APEX, idx(700.0), inc(pil_world_time(w))) * 180.0 /
         M_PI;
}

int main(int argc, char **argv) {
  static pr_t p;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  p.w = w;
  pil_world_add_overlay(w, overlay, &p);

  pil_plot *q = pil_add_plot(w, 0.64, 0.25, 0.99, 0.75, "\\theta_i [deg]",
                             "deviation [deg]", TSWEEP);
  pil_plot_add_signal(q, "", sig_deg, NULL);
  pil_plot_add_signal(q, "violet", sig_v, NULL);
  pil_plot_add_signal(q, "red", sig_r, NULL);
  pil_plot_set_xy(q, 0);
  pil_plots_record(w, 0.02, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -7.0, -2.4, 6.5, 3.4, 1300, 760, 0.0, 0.0, 0.62,
                      1.0, 14);
  int rc = example_main(w, &view, "Prism dispersion", argc, argv, "prism.png");

  pil_world_free(w);

  return rc;
}
