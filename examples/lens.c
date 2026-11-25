/*
 * lens.c - ray diagram of a thin converging lens (focal length f = 1)
 * Object (blue arrow) moves towards lens
 * 3 principal rays leave its tip: one parallel to axis (it is bent through far
 * focus), one through centre of lens (undeflected) and one through near focus
 * (it emerges parallel). Where they meet is image (red arrow): real and
 * inverted for do > f, at infinity for do = f, virtual, upright and magnified
 * (dashed back-extensions) for do < f. Plot: 1/d_i against 1/d_o, straight line
 * 1/d_i = 1/f - 1/d_o
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 3.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_optics.h>

#define F 1.0
#define H 0.55
#define TSWEEP 10.0
#define XMAX 7.0

typedef struct {
  const pil_world *w;
} le_t;

static double dobj(double t) {
  double s = t / TSWEEP;
  if (s > 1)
    s = 1;

  return 3.4 - 2.9 * s; // 3.4 down to 0.5
}

static void seg(cairo_t *cr, const pil_view *v, double x0, double y0, double x1,
                double y1) {
  pil_vec2 p = pil_view_to_screen(v, pil_v2(x0, y0));
  pil_vec2 q = pil_view_to_screen(v, pil_v2(x1, y1));
  cairo_move_to(cr, p.x, p.y);
  cairo_line_to(cr, q.x, q.y);
  cairo_stroke(cr);
}

static void arrow(cairo_t *cr, const pil_view *v, double x, double y,
                  double lw) {
  seg(cr, v, x, 0, x, y);
  double hd = y > 0 ? -0.1 : 0.1;
  seg(cr, v, x, y, x - 0.07, y + hd);
  seg(cr, v, x, y, x + 0.07, y + hd);
  (void)lw;
}

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const le_t *l = user;
  pil_style st = pil_style_textbook();
  double d = dobj(pil_world_time(l->w)), s[3];
  pil_lens_principal_slopes(F, d, H, s);
  double di = fabs(d - F) < 1e-9 ? 1e9 : pil_lens_image(F, d);
  double m = pil_lens_magnification(F, d);

  cairo_new_path(cr);
  cairo_set_line_width(cr, st.thin_w);
  cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
  seg(cr, v, -XMAX, 0, XMAX, 0); // optical axis
  cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
  cairo_set_line_width(cr, st.line_w);
  seg(cr, v, 0, -1.5, 0, 1.5); // lens (a line with arrowheads)
  seg(cr, v, 0, 1.5, -0.1, 1.35);
  seg(cr, v, 0, 1.5, 0.1, 1.35);
  seg(cr, v, 0, -1.5, -0.1, -1.35);
  seg(cr, v, 0, -1.5, 0.1, -1.35);
  for (int k = -1; k <= 1; k += 2) {
    pil_vec2 p = pil_view_to_screen(v, pil_v2(k * F, 0));
    cairo_arc(cr, p.x, p.y, 3.5, 0, 2 * M_PI);
    cairo_fill(cr);
    pil_draw_text(cr, (pil_vec2){p.x, p.y + 26}, "F", 0.5, &st);
  }

  // Principal rays: to lens, then after lens
  cairo_set_source_rgba(cr, 0.8, 0.1, 0.1, 0.9);
  cairo_set_line_width(cr, 1.5);
  for (int k = 0; k < 3; k++) {
    double yl = H + s[k] * d, s2 = pil_lens_ray_slope(F, yl, s[k]);
    seg(cr, v, -d, H, 0, yl);
    seg(cr, v, 0, yl, XMAX, yl + s2 * XMAX);
    if (di < 0) { // Virtual image: dashed back extension
      double dash[2] = {5, 5};
      cairo_set_dash(cr, dash, 2, 0);
      seg(cr, v, 0, yl, di, yl + s2 * di);
      cairo_set_dash(cr, NULL, 0, 0);
    }
  }

  // Object and image
  cairo_set_source_rgb(cr, 0.1, 0.3, 0.85);
  cairo_set_line_width(cr, 3.0);
  arrow(cr, v, -d, H, 3.0);
  if (di < 50) {
    cairo_set_source_rgb(cr, 0.8, 0.1, 0.1);
    arrow(cr, v, di, m * H, 3.0);
  }

  char buf[64];
  snprintf(buf, sizeof buf, "d_o = %.2f f", d / F);
  pil_draw_text(cr, pil_view_to_screen(v, pil_v2(-XMAX + 0.2, 1.55)), buf, 0.0,
                &st);
  if (di < 50) {
    snprintf(buf, sizeof buf, "d_i = %.2f f   m = %.2f", di / F, m);
  } else {
    snprintf(buf, sizeof buf, "image at infinity");
  }
  pil_draw_text(cr, pil_view_to_screen(v, pil_v2(-XMAX + 0.2, 1.2)), buf, 0.0,
                &st);
  cairo_new_path(cr);
}

static double sig_xo(const pil_world *w, void *u) {
  (void)u;

  return F / dobj(pil_world_time(w));
}

static double sig_xi(const pil_world *w, void *u) {
  (void)u;

  return 1.0 - F / dobj(pil_world_time(w)); // f / d_i
}

int main(int argc, char **argv) {
  static le_t l;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  l.w = w;
  pil_world_add_overlay(w, overlay, &l);

  pil_plot *p =
      pil_add_plot(w, 0.68, 0.28, 0.99, 0.72, "f / d_o", "f / d_i", TSWEEP);
  pil_plot_add_signal(p, "", sig_xo, NULL);
  pil_plot_add_signal(p, "f/d_i", sig_xi, NULL);
  pil_plot_set_xy(p, 0);
  pil_plots_record(w, 0.02, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -7.5, -2.0, 7.5, 2.0, 1300, 760, 0.0, 0.0, 0.66,
                      1.0, 14);
  int rc =
      example_main(w, &view, "Thin lens ray diagram", argc, argv, "lens.png");

  pil_world_free(w);

  return rc;
}
