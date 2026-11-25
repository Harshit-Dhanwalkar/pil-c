/*
 * snell.c - refraction, total internal reflection and Fresnel reflectance
 * A ray inside glass (n = 1.5) hits glass-air surface at an incidence angle
 * that grows from 0 to 85 degrees
 * Refraction obeys Snell's law n_1 sin(theta_i) = n_2 sin(theta_t); beyond
 * critical angle (41.8 degrees) ray is totally reflected
 * Line widths show how power is shared between reflected and transmitted rays
 * (unpolarised light)
 * Plot: Fresnel reflectance for two polarisations against incidence angle:
 * p-curve touches zero at Brewster's angle (33.7 degrees), both reach 1 at
 * critical angle
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 4.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_optics.h>

#define N1 1.5
#define N2 1.0
#define TSWEEP 10.0
#define MAXDEG 85.0
#define LEN 4.2

typedef struct {
  const pil_world *w;
} sn_t;

static double theta_i(double t) {
  double s = t / TSWEEP;

  return (s > 1 ? 1 : s) * MAXDEG * M_PI / 180.0;
}

static void ray(cairo_t *cr, const pil_view *v, double x0, double y0, double x1,
                double y1, double width, double a) {
  pil_vec2 p = pil_view_to_screen(v, pil_v2(x0, y0));
  pil_vec2 q = pil_view_to_screen(v, pil_v2(x1, y1));
  cairo_set_source_rgba(cr, 0.1, 0.1, 0.1, a);
  cairo_set_line_width(cr, width);
  cairo_move_to(cr, p.x, p.y);
  cairo_line_to(cr, q.x, q.y);
  cairo_stroke(cr);
}

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const sn_t *s = user;
  pil_style st = pil_style_textbook();
  double ti = theta_i(pil_world_time(s->w)), tt;
  int ok = pil_snell(N1, N2, ti, &tt);
  double r = 0.5 * (pil_fresnel_rs(N1, N2, ti) + pil_fresnel_rp(N1, N2, ti));

  // Glass below surface y = 0
  pil_vec2 a = pil_view_to_screen(v, pil_v2(-LEN - 0.5, 0));
  pil_vec2 b = pil_view_to_screen(v, pil_v2(LEN + 0.5, -LEN));
  cairo_new_path(cr);
  cairo_set_source_rgba(cr, 0.55, 0.75, 0.95, 0.35);
  cairo_rectangle(cr, a.x, a.y, b.x - a.x, b.y - a.y);
  cairo_fill(cr);
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_set_line_width(cr, st.line_w);
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, b.x, a.y);
  cairo_stroke(cr);

  // Normal (dashed)
  double dash[2] = {6, 5};
  pil_vec2 n0 = pil_view_to_screen(v, pil_v2(0, -LEN));
  pil_vec2 n1 = pil_view_to_screen(v, pil_v2(0, LEN));
  cairo_set_dash(cr, dash, 2, 0);
  cairo_set_source_rgb(cr, 0.4, 0.4, 0.4);
  cairo_set_line_width(cr, st.thin_w);
  cairo_move_to(cr, n0.x, n0.y);
  cairo_line_to(cr, n1.x, n1.y);
  cairo_stroke(cr);
  cairo_set_dash(cr, NULL, 0, 0);

  // Incident, reflected and transmitted rays
  double d[2] = {sin(ti), cos(ti)};
  double nrm[2] = {0, -1};
  double rd[2];
  double td[2]; // Incident travels up-right
  d[1] = cos(ti);
  pil_reflect(d, nrm, rd);
  ray(cr, v, -LEN * sin(ti), -LEN * cos(ti), 0, 0, 3.0, 1.0);
  ray(cr, v, 0, 0, LEN * rd[0], LEN * rd[1], 0.5 + 3.5 * r, 0.4 + 0.6 * r);
  if (pil_refract(d, nrm, N1, N2, td)) {
    ray(cr, v, 0, 0, LEN * td[0], LEN * td[1], 0.5 + 3.5 * (1 - r),
        0.4 + 0.6 * (1 - r));
  }
  (void)ok;
  (void)tt;

  char buf[64];
  snprintf(buf, sizeof buf, "\\theta_i = %.1f", ti * 180.0 / M_PI);
  pil_draw_text(cr, pil_view_to_screen(v, pil_v2(-LEN, 3.4)), buf, 0.0, &st);
  if (ok) {
    snprintf(buf, sizeof buf, "\\theta_t = %.1f", tt * 180.0 / M_PI);
  } else {
    snprintf(buf, sizeof buf, "total internal reflection");
  }
  pil_draw_text(cr, pil_view_to_screen(v, pil_v2(-LEN, 2.9)), buf, 0.0, &st);
  pil_draw_text(cr, pil_view_to_screen(v, pil_v2(LEN - 0.2, 0.35)),
                "air  n = 1", 1.0, &st);
  pil_draw_text(cr, pil_view_to_screen(v, pil_v2(LEN - 0.2, -0.65)),
                "glass  n = 1.5", 1.0, &st);
  cairo_new_path(cr);
}

static double sig_deg(const pil_world *w, void *u) {
  (void)u;

  return theta_i(pil_world_time(w)) * 180.0 / M_PI;
}

static double sig_rs(const pil_world *w, void *u) {
  (void)u;

  return pil_fresnel_rs(N1, N2, theta_i(pil_world_time(w)));
}

static double sig_rp(const pil_world *w, void *u) {
  (void)u;

  return pil_fresnel_rp(N1, N2, theta_i(pil_world_time(w)));
}

int main(int argc, char **argv) {
  static sn_t s;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  s.w = w;
  pil_world_add_overlay(w, overlay, &s);

  pil_plot *p = pil_add_plot(w, 0.62, 0.25, 0.99, 0.75, "\\theta_i [deg]",
                             "reflectance", TSWEEP);
  pil_plot_add_signal(p, "", sig_deg, NULL);
  pil_plot_add_signal(p, "R_s", sig_rs, NULL);
  pil_plot_add_signal(p, "R_p", sig_rp, NULL);
  pil_plot_set_xy(p, 0);
  pil_plots_record(w, 0.02, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -5.0, -4.6, 5.0, 4.6, 1300, 760, 0.0, 0.0, 0.6,
                      1.0, 14);
  int rc = example_main(w, &view, "Snell's law and Fresnel reflectance", argc,
                        argv, "snell.png");

  pil_world_free(w);

  return rc;
}
