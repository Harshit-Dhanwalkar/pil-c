/*
 * belts.c - charged particle trapped in a planetary dipole field
 * Dipole field B = B_0 R^3 / r^5 (-3 z x, -3 z y, r^2 - 3 z^2), northward at
 * equator like Earth's
 * A proton started on equator at r = 3 R with a pitch angle of 63 degrees does
 * three things at once, on three time scales: it gyrates about field line,
 * bounces between mirror points near 13 degrees of latitude and drifts slowly
 * around planet (westward for positive charge, from gradient and curvature of
 * field). Lower plots: height z and accumulated longitude
 * Trajectory by pil_particle3_step
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 70.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_lagrange.h>

#define B0 100.0
#define DT 0.004
#define TEND 150.0

typedef struct {
  double (*pos)[3];
  double *z;
  double *lon;
  int n;
  const pil_world *w;
  pil_proj3 pr;
} belts_t;

static void dipole(const double x[3], double t, double B[3], double E[3],
                   void *u) {
  (void)t;
  (void)u;
  double r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  double r5 = r2 * r2 * sqrt(r2);
  B[0] = -3.0 * B0 * x[2] * x[0] / r5;
  B[1] = -3.0 * B0 * x[2] * x[1] / r5;
  B[2] = B0 * (r2 - 3.0 * x[2] * x[2]) / r5;
  E[0] = E[1] = E[2] = 0.0;
}

static int idx(const belts_t *b) {
  int i = (int)(pil_world_time(b->w) / DT);

  return i < b->n ? i : b->n - 1;
}

static void curve(cairo_t *cr, const pil_view *v, const pil_proj3 *pr, double L,
                  double lon, double lat0, double lat1) {
  for (int i = 0; i <= 60; i++) {
    double lat = lat0 + (lat1 - lat0) * i / 60, r = L * cos(lat) * cos(lat);
    pil_vec2 q = pil_view_to_screen(
        v, pil_proj3_point(pr, r * cos(lat) * cos(lon), r * cos(lat) * sin(lon),
                           r * sin(lat)));
    if (i == 0) {
      cairo_move_to(cr, q.x, q.y);
    } else {
      cairo_line_to(cr, q.x, q.y);
    }
  }

  cairo_stroke(cr);
}

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const belts_t *b = user;
  pil_style st = pil_style_textbook();
  int cur = idx(b);

  // Planet: disc with equator and a meridian
  pil_vec2 O = pil_view_to_screen(v, pil_proj3_point(&b->pr, 0, 0, 0));
  double sc = pil_view_to_screen(v, pil_proj3_point(&b->pr, 0, 0, 1)).y;
  sc = fabs(sc - O.y) / cos(b->pr.pitch > 0 ? b->pr.pitch : 0.0);
  cairo_new_path(cr);
  cairo_arc(cr, O.x, O.y, 1.0 * sc * 0.97, 0, 2 * M_PI);
  cairo_close_path(cr);
  cairo_set_source_rgb(cr, 0.93, 0.93, 0.93);
  cairo_fill_preserve(cr);
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_set_line_width(cr, st.line_w);
  cairo_stroke(cr);
  cairo_new_path(cr);

  // Field lines of shells L = 2, 3 and 4 in four meridian planes
  cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
  cairo_set_line_width(cr, st.thin_w);
  double Ls[3] = {2.0, 3.0, 4.0};
  for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 3; j++) {
      double lmax = acos(sqrt(1.0 / Ls[j]));
      curve(cr, v, &b->pr, Ls[j], k * M_PI / 2, -lmax, lmax);
    }
  }

  // Wake: last 12 s, fading
  int span = (int)(12.0 / DT);
  int first = cur - span > 0 ? cur - span : 0;
  cairo_set_line_width(cr, 1.0);
  for (int i = first + 1; i <= cur; i += 2) {
    double age = (double)(cur - i) / span, g = 0.05 + 0.85 * age;
    pil_vec2 p0 = pil_view_to_screen(
        v, pil_proj3_point(&b->pr, b->pos[i - 1][0], b->pos[i - 1][1],
                           b->pos[i - 1][2]));
    pil_vec2 p1 = pil_view_to_screen(
        v, pil_proj3_point(&b->pr, b->pos[i][0], b->pos[i][1], b->pos[i][2]));
    cairo_set_source_rgb(cr, g, g, g);
    cairo_move_to(cr, p0.x, p0.y);
    cairo_line_to(cr, p1.x, p1.y);
    cairo_stroke(cr);
  }

  pil_vec2 q =
      pil_view_to_screen(v, pil_proj3_point(&b->pr, b->pos[cur][0],
                                            b->pos[cur][1], b->pos[cur][2]));
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_arc(cr, q.x, q.y, 5, 0, 2 * M_PI);
  cairo_fill(cr);
  pil_draw_text(cr, pil_v2(O.x, O.y - 1.0 * sc - 8), "N", 0.5, &st);
  cairo_new_path(cr);
}

static double sig_z(const pil_world *w, void *u) {
  (void)w;
  const belts_t *b = u;

  return b->z[idx(b)];
}

static double sig_lon(const pil_world *w, void *u) {
  (void)w;
  const belts_t *b = u;

  return b->lon[idx(b)] / (2.0 * M_PI);
}

int main(int argc, char **argv) {
  static belts_t b;
  b.n = (int)(TEND / DT);
  b.pos = malloc((size_t)b.n * sizeof *b.pos);
  b.z = malloc((size_t)b.n * sizeof *b.z);
  b.lon = malloc((size_t)b.n * sizeof *b.lon);

  // Guiding centre on equator at r = 3 R (proton starts one gyro-radius
  // inside), v_perp = 1, v_par = 0.5 along field
  pil_particle3 p = {{2.79, 0, 0}, {0, 1.0, 0.5}, 1.0};
  double prev = 0, acc = 0;
  for (int i = 0; i < b.n; i++) {
    for (int k = 0; k < 3; k++) {
      b.pos[i][k] = p.x[k];
    }

    double ph = atan2(p.x[1], p.x[0]);
    double d = ph - prev;
    if (i > 0) {
      if (d > M_PI) {
        d -= 2 * M_PI;
      }
      if (d < -M_PI) {
        d += 2 * M_PI;
      }
      acc += d;
    }

    prev = ph;
    b.z[i] = p.x[2];
    b.lon[i] = acc;
    pil_particle3_step(&p, dipole, NULL, i * DT, DT);
  }

  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  b.w = w;
  b.pr = pil_proj3_make(0, 0, 0.7, 0.45, 1.0);
  pil_world_add_overlay(w, overlay, &b);

  pil_plot *a =
      pil_add_plot(w, 0.62, 0.06, 0.99, 0.46, "t [s]", "z [R]", 100.0);
  pil_plot_add_signal(a, "z", sig_z, &b);
  pil_plot *l = pil_add_plot(w, 0.62, 0.56, 0.99, 0.97, "t [s]",
                             "longitude [turns]", 100.0);
  pil_plot_add_signal(l, "\\phi", sig_lon, &b);
  pil_plots_record(w, 0.1, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -4.3, -3.0, 4.3, 3.0, 1300, 760, 0.0, 0.0, 0.6,
                      1.0, 14);
  int rc = example_main(w, &view, "Particle trapped in a dipole field", argc,
                        argv, "belts.png");

  pil_world_free(w);
  free(b.pos);
  free(b.z);
  free(b.lon);

  return rc;
}
