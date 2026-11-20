/*
 * phase.c - damped pendulum and its phase portrait
 * A pendulum released at large angle swings with a period that grows with
 * amplitude and loses energy to drag
 * Right: theta' against theta, drawn with pil_plot_set_xy; dot follows spiral
 * in toward stable point at origin while angle against time is plotted below
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 7.0
#include "example_main.h"
#include <math.h>

typedef struct {
  pil_body *b;
  pil_vec2 pivot;
  double L;
} pend;

static double theta_of(const pend *p) {
  return atan2(p->b->pos.x - p->pivot.x, p->pivot.y - p->b->pos.y);
}

static double sig_theta(const pil_world *w, void *u) {
  (void)w;

  return theta_of(u);
}

static double sig_omega(const pil_world *w, void *u) {
  (void)w;
  const pend *p = u;
  double th = theta_of(p);

  return (p->b->vel.x * cos(th) + p->b->vel.y * sin(th)) / p->L;
}

int main(int argc, char **argv) {
  static pend p = {.pivot = {0, 2.0}, .L = 1.2};
  pil_world *w = pil_world_new();
  const double th0 = 2.6;

  pil_add_ground(w, (pil_vec2){-0.5, p.pivot.y}, (pil_vec2){0.5, p.pivot.y});
  p.b = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.09),
                                         .mass = 1.0,
                                         .pos = {p.pivot.x + p.L * sin(th0),
                                                 p.pivot.y - p.L * cos(th0)},
                                         .flags = PIL_HIDE_COM});
  pil_add_distance(w, pil_at(p.pivot.x, p.pivot.y), pil_on(p.b, 0, 0), p.L, 0);
  pil_add_drag(w, p.b, 0.18, 0.0);

  pil_plot *ph = pil_add_plot(w, 0.50, 0.04, 0.99, 0.60, "\xce\xb8 [rad]",
                              "\xce\xb8' [rad/s]", 14.0);
  pil_plot_add_signal(ph, "", sig_theta, &p); // x signal
  pil_plot_add_signal(ph, "", sig_omega, &p);
  pil_plot_set_xy(ph, 0);
  pil_plot *pt =
      pil_add_plot(w, 0.50, 0.64, 0.99, 0.98, "t [s]", "\xce\xb8 [rad]", 14.0);
  pil_plot_add_signal(pt, "\xce\xb8", sig_theta, &p);
  pil_plots_record(w, 0.02, 5e-4);

  pil_view view;
  pil_view_fit_region(&view, -1.6, -0.3, 1.6, 2.4, 1500, 800, 0.0, 0.0, 0.48,
                      1.0, 12);
  int rc = example_main(w, &view, "Damped pendulum: phase portrait", argc, argv,
                        "phase.png");

  pil_world_free(w);

  return rc;
}
