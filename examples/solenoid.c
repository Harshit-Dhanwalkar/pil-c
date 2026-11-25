/*
 * solenoid.c - magnetic field of a finite solenoid (Biot-Savart)
 * A solenoid of radius R and 4 turns per R is lengthened from L = 2 R to L = 12
 * R Field lines in a plane through axis are re-traced from summed Biot-Savart
 * contribution of every turn
 * Each field line carries a small arrowhead showing local direction of B; axial
 * line is highlighted in blue to show that B stays on axis and grows uniform as
 * solenoid lengthens
 * Axis labels (z, rho), pole labels (S left, N right - B exits +z end) and a
 * running length dimension
 * L = N R make connection between animation and plot below
 * Coil is drawn with dots on one side (current out of page) and crosses on
 * other (into page), following right-hand rule for B along +z inside
 * Swap two `side == +-1` branches if your convention has B reversed
 * Plot: field at centre divided by ideal-solenoid value 4 pi n climbs towards 1
 * with L / R, following L / sqrt(L^2 + 4 R^2)
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 8.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_biot.h>
#include <pil/pil_solenoid.h>
#include <stdio.h>

#define NPER 4.0
#define NSEG 36
#define TSPAN 10.0

typedef struct {
  const pil_world *w;
} s_t;

static double length(double t) {
  double s = t / TSPAN;

  return 2.0 + 10.0 * (s > 1 ? 1 : s);
}

typedef struct {
  double L;
  int n;
} sol_t;

static void field(double x, double z, double *bx, double *bz, void *user) {
  const sol_t *s = user;
  pil_solenoid(1.0, s->L, s->n, x, z, NSEG, bx, bz);
}

// Small filled arrowhead pointing from `from` to `to` (screen space)
static void arrow_head(cairo_t *cr, pil_vec2 from, pil_vec2 to, double size) {
  double dx = to.x - from.x, dy = to.y - from.y;
  double l = hypot(dx, dy);
  if (l < 1e-9) {
    return;
  }
  dx /= l;
  dy /= l;
  double nx = -dy, ny = dx;
  cairo_new_path(cr);
  cairo_move_to(cr, to.x, to.y);
  cairo_line_to(cr, to.x - size * dx + size * 0.45 * nx,
                to.y - size * dy + size * 0.45 * ny);
  cairo_line_to(cr, to.x - size * dx - size * 0.45 * nx,
                to.y - size * dy - size * 0.45 * ny);
  cairo_close_path(cr);
  cairo_fill(cr);
}

// Centered text label at screen position `p`
static void label(cairo_t *cr, pil_vec2 p, const char *txt, double size,
                  double r, double g, double b) {
  cairo_new_path(cr);
  cairo_select_font_face(cr, "sans-serif", CAIRO_FONT_SLANT_NORMAL,
                         CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, size);
  cairo_set_source_rgb(cr, r, g, b);
  cairo_text_extents_t e;
  cairo_text_extents(cr, txt, &e);
  cairo_move_to(cr, p.x - 0.5 * e.width - e.x_bearing,
                p.y + 0.5 * e.height - e.y_bearing);
  cairo_show_text(cr, txt);
}

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const s_t *g = user;
  pil_style st = pil_style_textbook();
  sol_t s = {length(pil_world_time(g->w)), 0};
  s.n = (int)lround(NPER * s.L);
  static double px[801], pz[801];
  double zend = 0.5 * s.L;

  // Field lines
  cairo_set_line_width(cr, st.thin_w * 1.3);
  double seeds[] = {0.0, 0.2, 0.4, 0.6, 1.5, 2.1};
  for (int i = 0; i < 6; i++) {
    for (int side = -1; side <= 1; side += 2) {
      for (int dir = -1; dir <= 1; dir += 2) {
        if (seeds[i] == 0.0 && side == 1) {
          continue;
        }
        int n = pil_field_line(field, &s, side * seeds[i], 0.0, dir, 0.05, 800,
                               2.8, 7.0, px, pz);
        int is_axis = (seeds[i] == 0.0);
        if (is_axis) {
          cairo_set_source_rgb(cr, 0.18, 0.38, 0.72);
        } else {
          cairo_set_source_rgb(cr, 0.25, 0.25, 0.25);
        }
        cairo_new_path(cr);
        for (int k = 0; k < n; k++) {
          pil_vec2 q = pil_view_to_screen(v, pil_v2(pz[k], px[k]));
          if (k == 0) {
            cairo_move_to(cr, q.x, q.y);
          } else {
            cairo_line_to(cr, q.x, q.y);
          }
        }
        cairo_stroke(cr);

        // Direction arrow at midpoint
        if (n > 30) {
          int a = n / 2;
          pil_vec2 f = pil_view_to_screen(v, pil_v2(pz[a - 1], px[a - 1]));
          pil_vec2 t = pil_view_to_screen(v, pil_v2(pz[a], px[a]));
          if (is_axis) {
            cairo_set_source_rgb(cr, 0.18, 0.38, 0.72);
          } else {
            cairo_set_source_rgb(cr, 0.25, 0.25, 0.25);
          }
          arrow_head(cr, f, t, 8.0);
        }
      }
    }
  }

  // Coil cross-sections
  // Dot  (side = +1, top of screen): current out of page
  // Cross(side = -1, bottom)       : current into page
  cairo_set_line_width(cr, 1.4);
  for (int k = 0; k < s.n; k++) {
    double z0 = -zend + (k + 0.5) * s.L / s.n;
    for (int side = -1; side <= 1; side += 2) {
      pil_vec2 q = pil_view_to_screen(v, pil_v2(z0, side * 1.0));
      cairo_set_source_rgb(cr, 0.82, 0.15, 0.12);
      cairo_new_path(cr);
      cairo_arc(cr, q.x, q.y, 2.8, 0, 2 * M_PI);
      if (side == 1) {
        cairo_fill(cr);
      } else {
        cairo_stroke(cr);
        cairo_new_path(cr);
        cairo_move_to(cr, q.x - 1.9, q.y - 1.9);
        cairo_line_to(cr, q.x + 1.9, q.y + 1.9);
        cairo_move_to(cr, q.x - 1.9, q.y + 1.9);
        cairo_line_to(cr, q.x + 1.9, q.y - 1.9);
        cairo_stroke(cr);
      }
    }
  }

  // Labels above and below coil
  {
    pil_vec2 qt = pil_view_to_screen(v, pil_v2(0.0, 1.42));
    label(cr, qt, "I", 13, 0.82, 0.15, 0.12);
    pil_vec2 qb = pil_view_to_screen(v, pil_v2(0.0, -1.42));
    label(cr, qb, "I", 13, 0.82, 0.15, 0.12);
  }

  // Axis labels
  {
    pil_vec2 qz = pil_view_to_screen(v, pil_v2(6.9, 0.3));
    label(cr, qz, "z", 16, 0.1, 0.1, 0.1);
    pil_vec2 qr = pil_view_to_screen(v, pil_v2(-6.9, 2.55));
    label(cr, qr, "ρ", 16, 0.1, 0.1, 0.1);
  }

  // Pole labels S (left), N (right)
  {
    pil_vec2 qs = pil_view_to_screen(v, pil_v2(-zend - 0.55, 1.65));
    label(cr, qs, "S", 18, 0.15, 0.30, 0.75);
    pil_vec2 qn = pil_view_to_screen(v, pil_v2(zend + 0.55, 1.65));
    label(cr, qn, "N", 18, 0.75, 0.15, 0.15);
  }

  // Length dimension above solenoid
  {
    double yt = 2.35;
    pil_vec2 a = pil_view_to_screen(v, pil_v2(-zend, yt));
    pil_vec2 b = pil_view_to_screen(v, pil_v2(zend, yt));
    cairo_set_source_rgb(cr, 0.4, 0.4, 0.4);
    cairo_set_line_width(cr, 1.0);
    cairo_new_path(cr);
    cairo_move_to(cr, a.x, a.y);
    cairo_line_to(cr, b.x, b.y);
    cairo_stroke(cr);
    for (int e = 0; e < 2; e++) {
      double zz = e ? zend : -zend;
      pil_vec2 t = pil_view_to_screen(v, pil_v2(zz, yt));
      cairo_new_path(cr);
      cairo_move_to(cr, t.x, t.y - 5);
      cairo_line_to(cr, t.x, t.y + 5);
      cairo_stroke(cr);
    }
    char buf[32];
    snprintf(buf, sizeof(buf), "L = %.0f R", s.L);
    pil_vec2 mid = {(a.x + b.x) / 2, a.y - 12};
    label(cr, mid, buf, 14, 0.25, 0.25, 0.25);
  }

  // Status: B_centre / (4 pi n)
  {
    char buf[80];
    double ratio =
        pil_solenoid_axis(1.0, s.L, NPER, 0.0) / pil_solenoid_ideal(NPER);
    snprintf(buf, sizeof(buf), "B_c / 4πn = %.3f", ratio);
    pil_vec2 q = pil_view_to_screen(v, pil_v2(-6.7, -2.5));
    label(cr, q, buf, 14, 0.25, 0.25, 0.25);
  }

  cairo_new_path(cr);
}

static double sig_l(const pil_world *w, void *u) {
  (void)u;

  return length(pil_world_time(w));
}

static double sig_b(const pil_world *w, void *u) {
  (void)u;
  double L = length(pil_world_time(w)), n = NPER;

  return pil_solenoid_axis(1.0, L, n, 0.0) / pil_solenoid_ideal(n);
}

int main(int argc, char **argv) {
  static s_t g;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  g.w = w;
  pil_world_add_overlay(w, overlay, &g);

  pil_plot *p = pil_add_plot(w, 0.68, 0.28, 0.99, 0.72, "L / R",
                             "B_{centre} / 4\\pi n", TSPAN);
  pil_plot_add_signal(p, "", sig_l, NULL);
  pil_plot_add_signal(p, "B(L)", sig_b, NULL);
  pil_plot_set_xy(p, 0);
  pil_plots_record(w, 0.05, 1e-3);

  pil_view view;
  pil_view_fit_region(&view, -7.0, -2.8, 7.0, 2.8, 1300, 760, 0.0, 0.0, 0.66,
                      1.0, 14);
  int rc =
      example_main(w, &view, "Finite solenoid", argc, argv, "solenoid.png");

  pil_world_free(w);

  return rc;
}
