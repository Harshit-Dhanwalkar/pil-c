/*
 * faraday.c - rectangular loop pulled through a magnetic field (B out of page)
 * Flux Phi = B h s changes only while loop enters (s growing) and leaves it
 * (s shrinking); emf = -d Phi / dt is then a pulse, of opposite sign for two,
 * and zero while whole loop is inside
 * Current I = emf / R feels force
 * B I h on edge in field, opposing motion (Lenz): loop is braked
 * only at entry and exit. Plot: Phi, emf and v
 * Braking force uses pil_add_custom_force, library's escape hatch
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 0.85
#include "example_main.h"
#include <math.h>

typedef struct {
  pil_body *loop;
  double B, h, l, W, R; // field, loop height/width, field width, resistance
} faraday;

static double clampd(double x, double a, double b) {
  return x < a ? a : (x > b ? b : x);
}

static int inside(const faraday *f, double x) { return x > 0.0 && x < f->W; }

// Loop edges at x -+ l/2: flux and its rate
static double flux(const faraday *f) {
  double xf = f->loop->pos.x + 0.5 * f->l, xr = f->loop->pos.x - 0.5 * f->l;

  return f->B * f->h * (clampd(xf, 0, f->W) - clampd(xr, 0, f->W));
}

static double flux_rate(const faraday *f) {
  double xf = f->loop->pos.x + 0.5 * f->l, xr = f->loop->pos.x - 0.5 * f->l;

  return f->B * f->h * f->loop->vel.x * (inside(f, xf) - inside(f, xr));
}

static void brake(pil_world *w, double t, void *user) {
  (void)w;
  (void)t;
  const faraday *f = user;
  // |I| = |emf| / R, force B |I| h on the edge in the field: by Lenz's law it
  // always opposes the motion, at entry and at exit alike
  double I = fabs(flux_rate(f)) / f->R;
  double dir = f->loop->vel.x >= 0 ? -1.0 : 1.0;
  pil_apply_force(f->loop, f->loop->pos, pil_v2(dir * f->B * f->h * I, 0));
}

static double sig_phi(const pil_world *w, void *user) {
  (void)w;

  return flux(user);
}

static double sig_emf(const pil_world *w, void *user) {
  (void)w;

  return -flux_rate(user);
}

static double sig_v(const pil_world *w, void *user) {
  (void)w;

  return ((const faraday *)user)->loop->vel.x;
}

int main(int argc, char **argv) {
  static faraday f = {.B = 1.0, .h = 1.0, .l = 1.6, .W = 3.0, .R = 2.5};
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);

  pil_add_field_marks(w, (pil_vec2){0.0, -0.9}, (pil_vec2){f.W, 0.9}, 0.4,
                      PIL_FIELD_B, 1, 0);

  const double a = 0.5 * f.l;
  const double b = 0.5 * f.h;
  const double t = 0.07;
  const double s = 0.004;
  pil_vec2 P[] = {{-a, -b},         {a, -b},        {a, b},
                  {-a, b},          {-a, s},        {-a + t, s},
                  {-a + t, b - t},  {a - t, b - t}, {a - t, -b + t},
                  {-a + t, -b + t}, {-a + t, -s},   {-a, -s}};
  const int n = 12;
  pil_vec2 c = pil_poly_centroid(P, n);
  for (int i = 0; i < n; i++) {
    P[i].x += -2.2 - c.x;
    P[i].y += -c.y;
  }

  c = pil_poly_centroid(P, n);
  f.loop = pil_add_body(w, &(pil_body_desc){.shape = pil_poly(P, n),
                                            .mass = 1.0,
                                            .pos = c,
                                            .vel = {3.0, 0},
                                            .fill = 0.75,
                                            .flags = PIL_HIDE_COM});
  pil_add_custom_force(w, brake, &f);
  pil_add_vector(w, pil_on(f.loop, 0, 0), PIL_VEC_VELOCITY, 30.0, 0, "v");
  pil_add_label(w, pil_at(0.0, 1.15), "B out of the page", 2, 0);
  pil_add_label(w, pil_at(f.W, 1.15), "loop: resistance R", 1, 0);

  pil_plot *pl = pil_add_plot(w, 0.08, 0.57, 0.98, 0.97, "t [s]",
                              "\xce\xa6 [Wb], emf [V], v [m/s]", 4.0);
  pil_plot_add_signal(pl, "\xce\xa6", sig_phi, &f);
  pil_plot_add_signal(pl, "emf", sig_emf, &f);
  pil_plot_add_signal(pl, "v", sig_v, &f);
  pil_plots_record(w, 0.02, 5e-4);

  pil_view view;
  pil_view_fit_region(&view, -1.6, -1.1, 4.6, 1.5, 1500, 900, 0.0, 0.0, 1.0,
                      0.6, 12);
  int rc =
      example_main(w, &view, "Faraday induction", argc, argv, "faraday.png");

  pil_world_free(w);

  return rc;
}
