/*
 * kapitza.c - Kapitza pendulum: vibrating pivot stabilises inverted pendulum
 * Both pendula start 23 degrees from upright position with their pivot shaking
 * vertically at w = 60 rad/s (shown in slow motion)
 * Left: amplitude a = 0.12 l, above critical a_c = sqrt(2 g l) / w = 0.074 l,
 * so rod stays up and wobbles slowly about vertical, in effective potential
 * V_eff = -(g/l) cos(theta) + a^2 w^2 sin^2(theta) / (4 l^2)
 * Right: a = 0.04 l is below it and pendulum falls over Plot shows angle of
 * each
 */

#define PIL_EXAMPLE_TIME_SCALE 0.1
#define PIL_EXAMPLE_IMAGE_TIME 1.6
#include "example_main.h"
#include <math.h>
#include <pil/pil_kapitza.h>
#include <pil/pil_ode.h>

#define TEND 3.0
#define H 5e-4

typedef struct {
  const pil_world *w;
} k_t;

static pil_kapitza par[2] = {{9.81, 1.0, 0.12, 60.0, 0.4},
                             {9.81, 1.0, 0.04, 60.0, 0.4}};
static double st[2][2], tnow = -1.0;

static void advance(const pil_world *w) {
  double t = pil_world_time(w);
  if (tnow < 0.0 || t < tnow) {
    for (int i = 0; i < 2; i++) {
      st[i][0] = 0.4;
      st[i][1] = 0.0;
    }
    tnow = 0.0;
  }

  while (t - tnow > 1e-12) {
    double h = t - tnow < H ? t - tnow : H;
    for (int i = 0; i < 2; i++) {
      pil_rk4_step(pil_kapitza_deriv, 2, tnow, st[i], h, &par[i]);
    }
    tnow += h;
  }
}

static void line(cairo_t *cr, const pil_view *v, double x0, double y0,
                 double x1, double y1) {
  pil_vec2 a = pil_view_to_screen(v, pil_v2(x0, y0));
  pil_vec2 b = pil_view_to_screen(v, pil_v2(x1, y1));
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, b.x, b.y);
  cairo_stroke(cr);
}

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const k_t *k = user;
  pil_style sty = pil_style_textbook();
  advance(k->w);
  double t = pil_world_time(k->w);
  for (int i = 0; i < 2; i++) {
    double cx = i == 0 ? -1.8 : 1.8;
    double yp = pil_kapitza_pivot(&par[i], t);
    double th = st[i][0];
    double bx = cx + par[i].l * sin(th);
    double by = yp + par[i].l * cos(th);
    cairo_new_path(cr);
    cairo_set_source_rgb(cr, 0.55, 0.55, 0.55);
    cairo_set_line_width(cr, sty.thin_w);
    line(cr, v, cx, -par[i].a, cx, par[i].a); // Range of pivot motion
    line(cr, v, cx - 0.4, -par[i].a, cx + 0.4, -par[i].a);
    line(cr, v, cx - 0.4, par[i].a, cx + 0.4, par[i].a);
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_set_line_width(cr, sty.line_w);
    line(cr, v, cx, yp, bx, by);
    pil_vec2 pv = pil_view_to_screen(v, pil_v2(cx, yp));
    pil_vec2 bv = pil_view_to_screen(v, pil_v2(bx, by));
    cairo_arc(cr, pv.x, pv.y, 4, 0, 2 * M_PI);
    cairo_fill(cr);
    cairo_set_source_rgb(cr, 0.85, 0.15, 0.1);
    cairo_arc(cr, bv.x, bv.y, 11, 0, 2 * M_PI);
    cairo_fill(cr);
    char buf[48];
    snprintf(buf, sizeof buf, "a = %.2f l  (%s)", par[i].a,
             par[i].a > pil_kapitza_critical_amp(par[i].g, par[i].l, par[i].w)
                 ? "a > a_c"
                 : "a < a_c");
    cairo_set_source_rgb(cr, 0, 0, 0);
    pil_draw_text(cr, pil_view_to_screen(v, pil_v2(cx, -1.45)), buf, 0.5, &sty);
  }

  cairo_new_path(cr);
}

static double sig_th(const pil_world *w, void *u) {
  advance(w);

  return st[*(int *)u][0] * 180.0 / M_PI;
}

int main(int argc, char **argv) {
  static k_t k;
  static int i0 = 0;
  static int i1 = 1;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  k.w = w;
  pil_world_add_overlay(w, overlay, &k);

  pil_plot *p =
      pil_add_plot(w, 0.62, 0.25, 0.99, 0.75, "t [s]", "\\theta [deg]", TEND);
  pil_plot_add_signal(p, "a > a_c", sig_th, &i0);
  pil_plot_add_signal(p, "a < a_c", sig_th, &i1);
  pil_plots_record(w, 0.004, 1.0);

  pil_view view;
  pil_view_fit_region(&view, -3.3, -1.7, 3.3, 1.5, 1300, 760, 0.0, 0.0, 0.6,
                      1.0, 14);
  int rc =
      example_main(w, &view, "Kapitza pendulum", argc, argv, "kapitza.png");

  pil_world_free(w);

  return rc;
}
