/*
 * rutherford.c - Rutherford scattering off a fixed repulsive nucleus
 * A beam of particles (speed v0, impact parameters b) is repelled by a charge
 * at origin with force K / r^2
 * Each path is a hyperbola; deflection follows tan(theta / 2) = kappa / b with
 * kappa = K / v0^2: small b gives large angles, head-on particle (b = 0 not
 * shown) is turned back at r_min = 2 kappa
 * Blue: simulated paths and particles (black dots); grey dashed: exit
 * asymptotes drawn with analytic angle, which paths follow
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 8.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_scatter.h>

#define NB 13
#define K_ 2.0
#define V0 2.0
#define X0 -9.0
#define DT 0.002
#define NMAX 20000

typedef struct {
  double (*path)[NMAX + 1][2];
  int n[NB];
  double b[NB];
  double theta[NB];
  const pil_world *w;
} ru_t;

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const ru_t *r = user;
  pil_style st = pil_style_textbook();
  int cur = (int)(pil_world_time(r->w) / DT);

  // Asymptotes of exit direction, dashed
  double dash[2] = {5, 5};
  cairo_new_path(cr);
  cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
  cairo_set_line_width(cr, st.thin_w);
  cairo_set_dash(cr, dash, 2, 0);
  for (int k = 0; k < NB; k++) {
    double th = pil_scatter_angle(K_, V0, fabs(r->b[k]));
    th = r->b[k] < 0 ? -th : th;
    double dx = cos(th);
    double dy = sin(th);
    double xe = r->path[k][r->n[k] - 1][0];
    double ye = r->path[k][r->n[k] - 1][1];
    double along = xe * dx + ye * dy;
    double px = xe - along * dx;
    double py = ye - along * dy; // Closest point of asymptote to nucleus
    pil_vec2 a = pil_view_to_screen(v, pil_v2(px, py));
    pil_vec2 b = pil_view_to_screen(v, pil_v2(px + 14 * dx, py + 14 * dy));
    cairo_move_to(cr, a.x, a.y);
    cairo_line_to(cr, b.x, b.y);
  }

  cairo_stroke(cr);
  cairo_set_dash(cr, NULL, 0, 0);

  // Paths so far and particles
  cairo_set_line_width(cr, 1.2);
  for (int k = 0; k < NB; k++) {
    int last = cur < r->n[k] ? cur : r->n[k] - 1;
    cairo_set_source_rgb(cr, 0.15, 0.3, 0.8);
    for (int i = 0; i < last; i += 5) {
      pil_vec2 q =
          pil_view_to_screen(v, pil_v2(r->path[k][i][0], r->path[k][i][1]));
      if (i == 0) {
        cairo_move_to(cr, q.x, q.y);
      } else {
        cairo_line_to(cr, q.x, q.y);
      }
    }

    pil_vec2 ql =
        pil_view_to_screen(v, pil_v2(r->path[k][last][0], r->path[k][last][1]));
    cairo_line_to(cr, ql.x, ql.y);

    cairo_stroke(cr);
    pil_vec2 q =
        pil_view_to_screen(v, pil_v2(r->path[k][last][0], r->path[k][last][1]));
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_arc(cr, q.x, q.y, 4, 0, 2 * M_PI);
    cairo_fill(cr);
  }

  // Nucleus
  pil_vec2 o = pil_view_to_screen(v, pil_v2(0, 0));
  cairo_set_source_rgb(cr, 0.8, 0.1, 0.1);
  cairo_arc(cr, o.x, o.y, 9, 0, 2 * M_PI);
  cairo_fill(cr);
  pil_draw_text(cr, (pil_vec2){o.x, o.y + 30}, "+", 0.5, &st);
  cairo_new_path(cr);
}

int main(int argc, char **argv) {
  static ru_t r;
  r.path = malloc(NB * sizeof *r.path);
  for (int k = 0; k < NB; k++) {
    int j = k - NB / 2; // -6 .. 6, never zero impact parameter
    r.b[k] = j == 0 ? 0.08 : 0.28 * j;
    r.n[k] = pil_scatter_trace(K_, V0, r.b[k], X0, DT, NMAX, r.path[k],
                               &r.theta[k], NULL);
  }

  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  r.w = w;
  pil_world_add_overlay(w, overlay, &r);

  pil_view view;
  pil_view_fit_region(&view, -10.0, -5.6, 10.0, 5.6, 1300, 760, 0.0, 0.0, 1.0,
                      1.0, 14);
  int rc = example_main(w, &view, "Rutherford scattering", argc, argv,
                        "rutherford.png");

  pil_world_free(w);
  free(r.path);

  return rc;
}
