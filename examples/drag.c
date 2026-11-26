/*
 * drag.c - projectile with quadratic air drag
 * Launch speed 10, k = 0.08 (a = -g y^ - k |v| v)
 * Launch angle sweeps from 10 to 80 degrees
 * Black: path with drag, grey dashed: same launch in vacuum (a parabola of
 * range v^2 sin(2 theta) / g)
 * Plot shows range of both against launch angle: drag shortens every throw,
 * steep ones most and best angle drops below 45 degrees
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 5.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_drag.h>

#define V0 10.0
#define TSW 10.0
#define MAXP 4000

static const pil_drag AIR = {9.81, 0.08};
static const VAC = {9.81, 0.0};

typedef struct {
  const pil_world *w;
} d_t;

static double ang(double t) {
  double s = t / TSW;

  return (10.0 + 70.0 * (s > 1 ? 1 : s)) * M_PI / 180.0;
}

static void path(cairo_t *cr, const pil_view *v, double (*pts)[2], int n) {
  for (int i = 0; i < n; i++) {
    pil_vec2 q = pil_view_to_screen(v, pil_v2(pts[i][0], pts[i][1]));
    if (i == 0) {
      cairo_move_to(cr, q.x, q.y);
    } else {
      cairo_line_to(cr, q.x, q.y);
    }
  }

  cairo_stroke(cr);
}

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const d_t *d = user;
  pil_style st = pil_style_textbook();
  double th = ang(pil_world_time(d->w));
  static double a[MAXP][2];
  static double b[MAXP][2];
  int na = pil_drag_trace(&AIR, 0.0, V0, th, 0.005, a, MAXP, NULL);
  int nb = pil_drag_trace(&VAC, 0.0, V0, th, 0.005, b, MAXP, NULL);

  cairo_new_path(cr);
  cairo_set_source_rgb(cr, 0.4, 0.4, 0.4);
  cairo_set_line_width(cr, st.thin_w);
  pil_vec2 g0 = pil_view_to_screen(v, pil_v2(-0.5, 0));
  pil_vec2 g1 = pil_view_to_screen(v, pil_v2(11.5, 0));
  cairo_move_to(cr, g0.x, g0.y);
  cairo_line_to(cr, g1.x, g1.y);
  cairo_stroke(cr);
  double dash[2] = {6, 5};
  cairo_set_source_rgb(cr, 0.55, 0.55, 0.55);
  cairo_set_dash(cr, dash, 2, 0);
  path(cr, v, b, nb);
  cairo_set_dash(cr, NULL, 0, 0);
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_set_line_width(cr, st.line_w);
  path(cr, v, a, na);
  pil_vec2 q = pil_view_to_screen(v, pil_v2(a[na - 1][0], 0));
  cairo_set_source_rgb(cr, 0.85, 0.15, 0.1);
  cairo_arc(cr, q.x, q.y, 5, 0, 2 * M_PI);
  cairo_fill(cr);
  char buf[48];
  snprintf(buf, sizeof buf, "\\theta = %.0f", th * 180.0 / M_PI);
  cairo_set_source_rgb(cr, 0, 0, 0);
  pil_draw_text(cr, pil_view_to_screen(v, pil_v2(0.0, 4.2)), buf, 0.0, &st);
  cairo_new_path(cr);
}

static double s_ang(const pil_world *w, void *u) {
  (void)u;
  return ang(pil_world_time(w)) * 180.0 / M_PI;
}

static double s_air(const pil_world *w, void *u) {
  (void)u;
  return pil_drag_range(&AIR, 0.0, V0, ang(pil_world_time(w)));
}

static double s_vac(const pil_world *w, void *u) {
  (void)u;
  return pil_drag_range(&VAC, 0.0, V0, ang(pil_world_time(w)));
}

int main(int argc, char **argv) {
  static d_t d;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  d.w = w;
  pil_world_add_overlay(w, overlay, &d);

  pil_plot *p =
      pil_add_plot(w, 0.68, 0.28, 0.99, 0.72, "\\theta [deg]", "range", TSW);
  pil_plot_add_signal(p, "", s_ang, NULL);
  pil_plot_add_signal(p, "drag", s_air, NULL);
  pil_plot_add_signal(p, "vacuum", s_vac, NULL);
  pil_plot_set_xy(p, 0);
  pil_plots_record(w, 0.05, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -0.8, -1.2, 11.8, 5.6, 1300, 760, 0.0, 0.0, 0.66,
                      1.0, 14);
  int rc = example_main(w, &view, "Projectile with air drag", argc, argv,
                        "drag.png");

  pil_world_free(w);

  return rc;
}
