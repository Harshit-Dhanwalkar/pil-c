/*
 * bottle.c - charged particle in a magnetic bottle
 * Axisymmetric field with B_z = B_0 (1 + z^2 / L^2) and, since div B = 0,
 * B_x = -B_0 x z / L^2, B_y = -B_0 y z / L^2: field lines pinch together near z
 * = 0 and crowd at ends
 * Particle spirals about a line and as field strengthens, its perpendicular
 * speed grows at expense of v_par until v_par = 0 and it is reflected (mirror
 * point)
 * Magnetic moment mu = m v_perp^2 / (2 B) is an adiabatic invariant: constant
 * to a few per cent along whole motion (lower plot)
 * Trajectory by pil_particle3_step
 * Everything Public is API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 30.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_lagrange.h>

#define B0 3.0
#define LB 3.0
#define DT 0.004
#define TEND 60.0

typedef struct {
  double (*pos)[3];
  double *z;
  double *mu;
  int n;
  const pil_world *w;
  pil_proj3 pr;
} bottle_t;

static void bfield(const double x[3], double t, double B[3], double E[3],
                   void *u) {
  (void)t;
  (void)u;
  B[0] = -B0 * x[0] * x[2] / (LB * LB);
  B[1] = -B0 * x[1] * x[2] / (LB * LB);
  B[2] = B0 * (1.0 + x[2] * x[2] / (LB * LB));
  E[0] = E[1] = E[2] = 0.0;
}

static double mu_of(const pil_particle3 *p) {
  double B[3];
  double E[3];
  bfield(p->x, 0, B, E, NULL);
  double bm = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
  double vpar = (p->v[0] * B[0] + p->v[1] * B[1] + p->v[2] * B[2]) / bm;
  double v2 = p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2];

  return 0.5 * (v2 - vpar * vpar) / bm;
}

static int idx(const bottle_t *b) {
  int i = (int)(pil_world_time(b->w) / DT);

  return i < b->n ? i : b->n - 1;
}

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const bottle_t *b = user;
  pil_style st = pil_style_textbook();
  int cur = idx(b);

  // Coils at z = +-2.4 and field lines: rho^2 B_z = const along a line
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_set_line_width(cr, st.line_w * 1.6);
  for (int s = -1; s <= 1; s += 2) {
    for (int i = 0; i <= 60; i++) {
      double a = 2 * M_PI * i / 60;
      pil_vec2 q = pil_view_to_screen(
          v, pil_proj3_point(&b->pr, 1.5 * cos(a), 1.5 * sin(a), s * 2.4));
      if (i == 0) {
        cairo_move_to(cr, q.x, q.y);
      } else {
        cairo_line_to(cr, q.x, q.y);
      }
    }

    cairo_stroke(cr);
  }

  cairo_set_source_rgb(cr, 0.55, 0.55, 0.55);
  cairo_set_line_width(cr, st.thin_w);
  for (int k = 0; k < 8; k++) {
    double a = 2 * M_PI * k / 8;
    for (int i = 0; i <= 40; i++) {
      double z = -2.4 + 4.8 * i / 40, rho = 1.0 / sqrt(1.0 + z * z / (LB * LB));
      pil_vec2 q = pil_view_to_screen(
          v, pil_proj3_point(&b->pr, rho * cos(a), rho * sin(a), z));
      if (i == 0) {
        cairo_move_to(cr, q.x, q.y);
      } else {
        cairo_line_to(cr, q.x, q.y);
      }
    }

    cairo_stroke(cr);
  }

  // Wake: last 12 s, fading
  int first = cur - (int)(12.0 / DT) > 0 ? cur - (int)(12.0 / DT) : 0;
  cairo_set_line_width(cr, 1.0);
  for (int i = first + 1; i <= cur; i += 2) {
    double age = (double)(cur - i) / (12.0 / DT);
    double g = 0.05 + 0.85 * age;
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
  pil_vec2 top = pil_view_to_screen(v, pil_proj3_point(&b->pr, 0, 0, 2.9));
  pil_draw_text(cr, pil_v2(top.x, top.y - 4), "B", 0.5, &st);
}

static double sig_z(const pil_world *w, void *u) {
  (void)w;
  const bottle_t *b = u;

  return b->z[idx(b)];
}

static double sig_mu(const pil_world *w, void *u) {
  (void)w;
  const bottle_t *b = u;

  return b->mu[idx(b)] / b->mu[0];
}

int main(int argc, char **argv) {
  static bottle_t b;
  b.n = (int)(TEND / DT);
  b.pos = malloc((size_t)b.n * sizeof *b.pos);
  b.z = malloc((size_t)b.n * sizeof *b.z);
  b.mu = malloc((size_t)b.n * sizeof *b.mu);

  pil_particle3 p = {
      {0.33, 0, 0}, {0, 1.0, 0.35}, 1.0}; // v_perp = 1, v_par = 0.35
  for (int i = 0; i < b.n; i++) {
    for (int k = 0; k < 3; k++) {
      b.pos[i][k] = p.x[k];
    }

    b.z[i] = p.x[2];
    b.mu[i] = mu_of(&p);
    pil_particle3_step(&p, bfield, NULL, i * DT, DT);
  }

  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  b.w = w;
  b.pr = pil_proj3_make(0, 0, 0.5, 0.25, 1.0);
  pil_world_add_overlay(w, overlay, &b);

  pil_plot *a = pil_add_plot(w, 0.62, 0.06, 0.99, 0.46, "t [s]", "z", 40.0);
  pil_plot_add_signal(a, "z", sig_z, &b);
  pil_plot *m =
      pil_add_plot(w, 0.62, 0.56, 0.99, 0.97, "t [s]", "\\mu / \\mu_0", 40.0);
  pil_plot_add_signal(m, "\\mu", sig_mu, &b);
  pil_plots_record(w, 0.05, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -1.8, -3.1, 1.8, 3.1, 1300, 760, 0.0, 0.0, 0.6,
                      1.0, 14);
  int rc = example_main(w, &view, "Magnetic bottle", argc, argv, "bottle.png");

  pil_world_free(w);
  free(b.pos);
  free(b.z);
  free(b.mu);

  return rc;
}
