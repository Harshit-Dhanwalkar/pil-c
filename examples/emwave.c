/*
 * emwave.c - a plane electromagnetic wave and its Poynting vector
 * A linearly polarised wave travels along x (c = 2 here). E (red) oscillates
 * along z, B (blue) along y, in phase and perpendicular to each other and to
 * direction of travel; B = E / c
 * Energy flux S = E x B / mu0 (black)
 * always points along +x and is never negative, pulsing at twice wave
 * frequency
 * Plot: energy density u(t) and S(t) at a fixed place, with two equal halves of
 * u (electric, magnetic) and mean intensity
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 2.2
#include "example_main.h"
#include <math.h>
#include <pil/pil_emwave.h>
#include <pil/pil_lagrange.h>

#define NS 49
#define XMAX 6.0

typedef struct {
  pil_emwave wv;
  pil_proj3 pr;
  const pil_world *w;
} em_t;

// arrow in the current source colour: shaft plus a small two-stroke head
static void arrow(cairo_t *cr, pil_vec2 a, pil_vec2 b, double head, double lw) {
  double dx = b.x - a.x, dy = b.y - a.y, l = hypot(dx, dy);
  if (l < 1.0) {
    return;
  }
  dx /= l;
  dy /= l;
  cairo_set_line_width(cr, lw);
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, b.x, b.y);
  cairo_move_to(cr, b.x - head * (dx + 0.45 * dy),
                b.y - head * (dy - 0.45 * dx));
  cairo_line_to(cr, b.x, b.y);
  cairo_line_to(cr, b.x - head * (dx - 0.45 * dy),
                b.y - head * (dy + 0.45 * dx));
  cairo_stroke(cr);
}

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const em_t *m = user;
  pil_style st = pil_style_textbook();
  double t = pil_world_time(m->w);

  // Propagation axis
  cairo_new_path(cr);
  cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
  cairo_set_line_width(cr, st.thin_w);
  pil_vec2 a = pil_view_to_screen(v, pil_proj3_point(&m->pr, -XMAX, 0, 0));
  pil_vec2 b = pil_view_to_screen(v, pil_proj3_point(&m->pr, XMAX + 0.6, 0, 0));
  cairo_move_to(cr, a.x, a.y);
  cairo_line_to(cr, b.x, b.y);
  cairo_stroke(cr);

  // E envelope (in the x-z plane) and B envelope (in the x-y plane)
  for (int k = 0; k < 2; k++) {
    if (k == 0) {
      cairo_set_source_rgba(cr, 0.8, 0.1, 0.1, 0.9);
    } else {
      cairo_set_source_rgba(cr, 0.1, 0.3, 0.85, 0.9);
    }
    cairo_set_line_width(cr, 1.4);
    for (int i = 0; i <= 240; i++) {
      double x = -XMAX + 2 * XMAX * i / 240.0;
      double E[3];
      double B[3];
      pil_emwave_field(&m->wv, x, t, E, B);
      double val = k == 0 ? E[1] : B[2] * m->wv.c;
      pil_vec2 q =
          pil_view_to_screen(v, k == 0 ? pil_proj3_point(&m->pr, x, 0, val)
                                       : pil_proj3_point(&m->pr, x, val, 0));
      if (i == 0) {
        cairo_move_to(cr, q.x, q.y);
      } else {
        cairo_line_to(cr, q.x, q.y);
      }
    }

    cairo_stroke(cr);
  }

  // Field arrows at sample points, drawn far to near
  for (int i = 0; i < NS; i++) {
    double x = -XMAX + 2 * XMAX * i / (NS - 1.0);
    double E[3];
    double B[3];
    double S[3];
    pil_emwave_field(&m->wv, x, t, E, B);
    pil_emwave_poynting(&m->wv, x, t, S);
    pil_vec2 o = pil_view_to_screen(v, pil_proj3_point(&m->pr, x, 0, 0));
    pil_vec2 pe = pil_view_to_screen(v, pil_proj3_point(&m->pr, x, 0, E[1]));
    pil_vec2 pb =
        pil_view_to_screen(v, pil_proj3_point(&m->pr, x, B[2] * m->wv.c, 0));
    pil_vec2 ps = pil_view_to_screen(
        v, pil_proj3_point(
               &m->pr, x + 0.8 * S[0] / (m->wv.c * m->wv.e0 * m->wv.e0), 0, 0));
    cairo_set_source_rgb(cr, 0.8, 0.1, 0.1);
    arrow(cr, o, pe, 6, 1.3);
    cairo_set_source_rgb(cr, 0.1, 0.3, 0.85);
    arrow(cr, o, pb, 6, 1.3);
    if (i % 3 == 0) { // Energy flux, drawn below axis
      cairo_set_source_rgb(cr, 0, 0, 0);
      pil_vec2 s0 = {o.x, o.y + 30}, s1 = {ps.x, o.y + 30};
      arrow(cr, s0, s1, 7, 2.4);
    }
  }

  // Key
  static const char *names[3] = {"E", "c B", "S (energy flux)"};
  double cols[3][3] = {{0.8, 0.1, 0.1}, {0.1, 0.3, 0.85}, {0, 0, 0}};
  for (int k = 0; k < 3; k++) {
    pil_vec2 p = pil_view_to_screen(v, pil_v2(-5.7, -2.2 - 0.45 * k));
    cairo_set_source_rgb(cr, cols[k][0], cols[k][1], cols[k][2]);
    cairo_set_line_width(cr, k == 2 ? 2.4 : 1.4);
    cairo_move_to(cr, p.x, p.y);
    cairo_line_to(cr, p.x + 26, p.y);
    cairo_stroke(cr);
    pil_draw_text(cr, (pil_vec2){p.x + 34, p.y + 5}, names[k], 0.0, &st);
  }

  pil_draw_text(
      cr, pil_view_to_screen(v, pil_proj3_point(&m->pr, XMAX + 0.8, 0, 0)), "x",
      0.0, &st);
  cairo_new_path(cr);
}

static double sig_u(const pil_world *w, void *u) {
  const em_t *m = u;

  return pil_emwave_energy_density(&m->wv, 0.0, pil_world_time(w));
}

static double sig_s(const pil_world *w, void *u) {
  const em_t *m = u;
  double S[3];
  pil_emwave_poynting(&m->wv, 0.0, pil_world_time(w), S);

  return S[0];
}

static double sig_i(const pil_world *w, void *u) {
  (void)w;
  const em_t *m = u;

  return pil_emwave_intensity(&m->wv);
}

int main(int argc, char **argv) {
  static em_t m;
  m.wv = (pil_emwave){2.0, 1.0, 2.0 * M_PI / 4.0, 0.0, 0.0};
  m.pr = pil_proj3_make(0, 0, 0.8, 0.25, 1.0);
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  m.w = w;
  pil_world_add_overlay(w, overlay, &m);

  pil_plot *p = pil_add_plot(w, 0.62, 0.20, 0.99, 0.80, "t", "u,  S", 8.0);
  pil_plot_add_signal(p, "S", sig_s, &m);
  pil_plot_add_signal(p, "u", sig_u, &m);
  pil_plot_add_signal(p, "<S>", sig_i, &m);
  pil_plots_record(w, 0.02, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -6.0, -3.2, 6.0, 3.2, 1300, 760, 0.0, 0.0, 0.6,
                      1.0, 14);
  int rc = example_main(w, &view, "Electromagnetic wave and Poynting vector",
                        argc, argv, "emwave.png");

  pil_world_free(w);

  return rc;
}
