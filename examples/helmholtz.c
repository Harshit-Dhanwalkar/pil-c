/*
 * helmholtz.c - magnetic field of two coaxial current loops (Biot-Savart)
 * Two loops of radius R carry same current
 * Separation d is slowly increased from 0.5 R to 1.5 R and field lines in plane
 * through axis are re-traced by summing Biot-Savart contribution of 200
 * segments of each loop. At d = R (Helmholtz pair) on-axis field is flat to
 * fourth order, so lines between coils run parallel
 * Plot: curvature of on-axis field:
 *  [B(h) + B(-h) - 2 B(0)] / (B(0) h^2) with h = 0.15 R, against d / R
 * It is negative for d < R (field peaked in middle), positive for d > R (dip in
 * middle) and crosses zero at d = R
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 5.0
#include "example_main.h"
#include <math.h>

#define NSEG 200
#define TSPAN 10.0

static double sep(double t) {
  double s = t / TSPAN;
  if (s > 1) {
    s = 1;
  }

  return 0.5 + s;
}

// B at (x, 0, z) from loops at z = +-d/2 (unit current, mu0/4pi = 1)
static void field(double d, double x, double z, double *bx, double *bz) {
  *bx = *bz = 0;
  for (int k = 0; k < 2; k++) {
    double z0 = k ? d / 2 : -d / 2;
    for (int i = 0; i < NSEG; i++) {
      double ph = 2 * M_PI * (i + 0.5) / NSEG, dph = 2 * M_PI / NSEG;
      double px = cos(ph);
      double py = sin(ph);
      double dlx = -py * dph;
      double dly = px * dph;
      double rx = x - px, ry = -py, rz = z - z0;
      double r3 = pow(rx * rx + ry * ry + rz * rz, 1.5);
      *bx += (dly * rz) / r3; // (dl x r)_x = dl_y r_z - dl_z r_y
      *bz += (dlx * ry - dly * rx) / r3;
    }
  }
}

typedef struct {
  const pil_world *w;
} hh_t;

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const hh_t *h = user;
  pil_style st = pil_style_textbook();
  double d = sep(pil_world_time(h->w));

  // field lines from seeds on the mid-plane, traced both ways
  cairo_new_path(cr);
  cairo_set_source_rgb(cr, 0.25, 0.25, 0.25);
  cairo_set_line_width(cr, st.thin_w * 1.3);
  double seeds[] = {0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.15, 1.4, 1.8};
  for (int s = 0; s < 10; s++) {
    for (int sgn = -1; sgn <= 1; sgn += 2) {
      for (int side = -1; side <= 1; side += 2) {
        if (seeds[s] == 0.0 && side == 1) {
          continue;
        }
        double x = side * seeds[s];
        double z = 0;
        pil_vec2 q = pil_view_to_screen(v, pil_v2(z, x));
        cairo_move_to(cr, q.x, q.y);
        for (int n = 0; n < 700; n++) {
          double bx;
          double bz;
          field(d, x, z, &bx, &bz);
          double bm = hypot(bx, bz);
          if (bm < 1e-9) {
            break;
          }
          // Midpoint step along the line
          double hs = 0.03 * sgn;
          double mx = x + 0.5 * hs * bx / bm;
          double mz = z + 0.5 * hs * bz / bm;
          field(d, mx, mz, &bx, &bz);
          bm = hypot(bx, bz);
          x += hs * bx / bm;
          z += hs * bz / bm;
          if (fabs(x) > 2.6 || fabs(z) > 3.4) {
            break;
          }
          q = pil_view_to_screen(v, pil_v2(z, x));
          cairo_line_to(cr, q.x, q.y);
          if (n > 20 && hypot(x - side * seeds[s], z) < 0.03) {
            break;
          }
        }
        cairo_stroke(cr);
      }
    }
  }

  // Loops seen edge-on: current out of the page at x = +R, into it at x = -R
  for (int k = 0; k < 2; k++) {
    double z0 = k ? d / 2 : -d / 2;
    for (int sd = -1; sd <= 1; sd += 2) {
      pil_vec2 q = pil_view_to_screen(v, pil_v2(z0, sd * 1.0));
      cairo_set_source_rgb(cr, 1, 1, 1);
      cairo_arc(cr, q.x, q.y, 9, 0, 2 * M_PI);
      cairo_fill_preserve(cr);
      cairo_set_source_rgb(cr, 0, 0, 0);
      cairo_set_line_width(cr, st.line_w);
      cairo_stroke(cr);
      if (sd > 0) {
        cairo_arc(cr, q.x, q.y, 2.5, 0, 2 * M_PI);
        cairo_fill(cr);
      } else {
        cairo_move_to(cr, q.x - 5, q.y - 5);
        cairo_line_to(cr, q.x + 5, q.y + 5);
        cairo_move_to(cr, q.x + 5, q.y - 5);
        cairo_line_to(cr, q.x - 5, q.y + 5);
        cairo_stroke(cr);
      }
    }
  }

  char buf[32];
  snprintf(buf, sizeof buf, "d = %.2f R", d);
  pil_draw_text(cr, pil_view_to_screen(v, pil_v2(0, -2.4)), buf, 0.5, &st);
  cairo_new_path(cr);
}

static double sig_dev(const pil_world *w, void *u) {
  (void)u;
  double d = sep(pil_world_time(w));
  double bx;
  double b0;
  double b1;
  double b2;
  double h = 0.15;
  field(d, 0, 0, &bx, &b0);
  field(d, 0, h, &bx, &b1);
  field(d, 0, -h, &bx, &b2);

  return (b1 + b2 - 2 * b0) / (b0 * h * h);
}

static double sig_d(const pil_world *w, void *u) {
  (void)u;

  return sep(pil_world_time(w));
}

int main(int argc, char **argv) {
  static hh_t h;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  h.w = w;
  pil_world_add_overlay(w, overlay, &h);

  pil_plot *p = pil_add_plot(w, 0.66, 0.30, 0.99, 0.70, "d / R",
                             "axis curvature of B_z [1/R^2]", TSPAN);
  pil_plot_add_signal(p, "", sig_d, NULL);
  pil_plot_add_signal(p, "curvature", sig_dev, NULL);
  pil_plot_set_xy(p, 0);
  pil_plots_record(w, 0.05, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -3.4, -2.6, 3.4, 2.6, 1300, 760, 0.0, 0.0, 0.62,
                      1.0, 14);
  int rc = example_main(w, &view, "Helmholtz coils: Biot-Savart field lines",
                        argc, argv, "helmholtz.png");
  pil_world_free(w);

  return rc;
}
