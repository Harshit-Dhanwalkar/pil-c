/*
 * resonance.c - driven damped oscillator, frequency swept slowly through
 * resonance Block m on spring k with dashpot c, pushed by F0 cos(phi(t)) whose
 * frequency rises linearly
 * Envelope of x(t) grows as drive nears natural frequency f0 = sqrt(k / m) / 2
 * pi and then dies away: plotted against drive frequency it is resonance curve,
 * peak height F0 / (c omega0) (quality factor Q = sqrt(k m) / c), width about
 * f0 / Q Drive uses pil_add_custom_force
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 10.5
#include "example_main.h"
#include <math.h>

typedef struct {
  pil_body *b;
  double F0, fa, beta; // drive: force amplitude, start frequency, df/dt
} drive;

static drive D;

static double drive_f(double t) { return D.fa + D.beta * t; }

static void push(pil_world *w, double t, void *user) {
  (void)w;
  drive *d = user;
  double phase = 2.0 * M_PI * (d->fa * t + 0.5 * d->beta * t * t);
  pil_apply_force(d->b, d->b->pos, pil_v2(d->F0 * cos(phase), 0));
}

static double sig_x(const pil_world *w, void *u) {
  (void)w;

  return ((const drive *)u)->b->pos.x - 1.5;
}

static double sig_f(const pil_world *w, void *u) {
  (void)u;

  return drive_f(pil_world_time(w));
}

// Amplitude of x from x and v at instantaneous drive frequency
static double sig_env(const pil_world *w, void *u) {
  const drive *d = u;
  double om = 2.0 * M_PI * drive_f(pil_world_time(w));
  double x = d->b->pos.x - 1.5, v = d->b->vel.x;

  return sqrt(x * x + v * v / (om * om));
}

int main(int argc, char **argv) {
  const double m = 1.0, k = 100.0, c = 2.0; // f0 = 1.59 Hz, Q = 5
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);

  D.F0 = 1.0;
  D.fa = 0.9;
  D.beta = 0.07; // 0.9 .. 2.3 Hz in 20 s
  D.b = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.5, 0.4),
                                         .mass = m,
                                         .pos = {1.5, 0.0},
                                         .fill = 0.85,
                                         .flags = PIL_HIDE_COM,
                                         .label = "m"});
  pil_add_prismatic(w, pil_on(D.b, 0, 0), pil_at(1.5, 0), (pil_vec2){1, 0});
  pil_add_ground(w, (pil_vec2){0, 0.5}, (pil_vec2){0, -0.7});
  pil_add_ground(w, (pil_vec2){-0.1, -0.2}, (pil_vec2){3.0, -0.2});
  pil_add_spring(w, pil_at(0, 0.1), pil_on(D.b, -0.25, 0.1), k, 1.25);
  pil_add_damper(w, pil_at(0, -0.1), pil_on(D.b, -0.25, -0.1), c);
  pil_add_custom_force(w, push, &D);
  pil_add_label(w, pil_at(0.6, 0.45), "k", 0, 0);
  pil_add_label(w, pil_at(0.6, -0.45), "c", 0, 0);

  pil_plot *px =
      pil_add_plot(w, 0.04, 0.34, 0.50, 0.98, "t [s]", "x [m]", 20.0);
  pil_plot_add_signal(px, "x", sig_x, &D);
  pil_plot *pr =
      pil_add_plot(w, 0.55, 0.34, 0.99, 0.98, "f [Hz]", "amplitude [m]", 20.0);
  pil_plot_add_signal(pr, "", sig_f, &D); // x axis
  pil_plot_add_signal(pr, "", sig_env, &D);
  pil_plot_set_xy(pr, 0);
  pil_plots_record(w, 0.02, 5e-4);

  pil_view view;
  pil_view_fit_region(&view, -0.3, -0.8, 3.1, 0.8, 1500, 900, 0.0, 0.0, 1.0,
                      0.28, 12);
  int rc =
      example_main(w, &view, "Resonance sweep", argc, argv, "resonance.png");

  pil_world_free(w);

  return rc;
}
