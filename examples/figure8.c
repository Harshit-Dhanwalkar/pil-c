/*
 * figure8.c - three bodies on a figure-eight orbit (Chenciner & Montgomery)
 * Three equal masses chase each other along one figure-eight curve, a periodic
 * solution of Newton's three-body problem found numerically by Moore (1993);
 * motion repeats after T = 6.3259 and every separation oscillates with period
 * T / 3
 * Plots: two of pairwise separations and energy error of leapfrog integrator,
 * which stays bounded
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 12.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_nbody.h>

#define DT 0.0005
#define TEND 40.0
#define TRAIL 3.2

typedef struct {
  double (*x)[3][2];
  double *r12;
  double *r23;
  double *de;
  int n;
  const pil_world *w;
} fig_t;

static int idx(const fig_t *f) {
  int i = (int)(pil_world_time(f->w) / DT);

  return i < f->n ? i : f->n - 1;
}

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const fig_t *f = user;
  pil_style st = pil_style_textbook();
  int cur = idx(f), span = (int)(TRAIL / DT);
  int first = cur - span > 0 ? cur - span : 0;
  double col[3][3] = {{0.8, 0.15, 0.15}, {0.15, 0.4, 0.85}, {0.1, 0.6, 0.25}};

  cairo_new_path(cr);
  cairo_set_line_width(cr, 1.6);
  for (int b = 0; b < 3; b++) {
    for (int i = first + 1; i <= cur; i += 4) {
      double a = 0.1 + 0.9 * (double)(i - first) / (cur - first + 1);
      pil_vec2 p0 =
          pil_view_to_screen(v, pil_v2(f->x[i - 1][b][0], f->x[i - 1][b][1]));
      pil_vec2 p1 = pil_view_to_screen(v, pil_v2(f->x[i][b][0], f->x[i][b][1]));
      cairo_set_source_rgba(cr, col[b][0], col[b][1], col[b][2], a);
      cairo_move_to(cr, p0.x, p0.y);
      cairo_line_to(cr, p1.x, p1.y);
      cairo_stroke(cr);
    }
  }

  for (int b = 0; b < 3; b++) {
    pil_vec2 q =
        pil_view_to_screen(v, pil_v2(f->x[cur][b][0], f->x[cur][b][1]));
    cairo_set_source_rgb(cr, col[b][0], col[b][1], col[b][2]);
    cairo_arc(cr, q.x, q.y, 9, 0, 2 * M_PI);
    cairo_fill_preserve(cr);
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_set_line_width(cr, st.thin_w);
    cairo_stroke(cr);
  }

  cairo_new_path(cr);
}

static double sig_r12(const pil_world *w, void *u) {
  (void)w;
  const fig_t *f = u;

  return f->r12[idx(f)];
}

static double sig_r23(const pil_world *w, void *u) {
  (void)w;
  const fig_t *f = u;

  return f->r23[idx(f)];
}

static double sig_de(const pil_world *w, void *u) {
  (void)w;
  const fig_t *f = u;

  return f->de[idx(f)];
}

int main(int argc, char **argv) {
  static fig_t f;
  f.n = (int)(TEND / DT);
  f.x = malloc((size_t)f.n * sizeof *f.x);
  f.r12 = malloc((size_t)f.n * sizeof *f.r12);
  f.r23 = malloc((size_t)f.n * sizeof *f.r23);
  f.de = malloc((size_t)f.n * sizeof *f.de);

  pil_nbody s;
  pil_nbody_figure_eight(&s);
  double e0 = pil_nbody_energy(&s);
  for (int i = 0; i < f.n; i++) {
    for (int b = 0; b < 3; b++) {
      f.x[i][b][0] = s.x[b][0];
      f.x[i][b][1] = s.x[b][1];
    }

    f.r12[i] = hypot(s.x[0][0] - s.x[1][0], s.x[0][1] - s.x[1][1]);
    f.r23[i] = hypot(s.x[1][0] - s.x[2][0], s.x[1][1] - s.x[2][1]);
    f.de[i] = pil_nbody_energy(&s) - e0;
    pil_nbody_step(&s, DT);
  }

  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  f.w = w;
  pil_world_add_overlay(w, overlay, &f);

  pil_plot *p =
      pil_add_plot(w, 0.62, 0.06, 0.99, 0.46, "t", "separation", TEND);
  pil_plot_add_signal(p, "r_{12}", sig_r12, &f);
  pil_plot_add_signal(p, "r_{23}", sig_r23, &f);
  pil_plot *q = pil_add_plot(w, 0.62, 0.56, 0.99, 0.97, "t", "\\Delta E", TEND);
  pil_plot_add_signal(q, "\\Delta E", sig_de, &f);
  pil_plots_record(w, 0.05, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -1.5, -0.9, 1.5, 0.9, 1300, 760, 0.0, 0.0, 0.6,
                      1.0, 14);
  int rc = example_main(w, &view, "Three bodies on a figure-eight orbit", argc,
                        argv, "figure8.png");

  pil_world_free(w);
  free(f.x);
  free(f.r12);
  free(f.r23);
  free(f.de);

  return rc;
}
