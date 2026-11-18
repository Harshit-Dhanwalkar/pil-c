/*
 * kepler.c - planet on an ellipse: Kepler's laws from an inverse-square force
 * A planet orbits a fixed star (GM = 1) with eccentricity 0.6
 * Orbit is traced, velocity and acceleration are drawn live and radius vector
 * sweeps a shaded wedge
 * Two fixed wedges of equal time (at perihelion and aphelion) have equal area;
 * plot shows distance, speed and areal velocity, which stays constant
 * Everything is public-API
 */

#define PIL_EXAMPLE_IMAGE_TIME 4.4
#define PIL_EXAMPLE_TIME_SCALE 0.70          // 4.4 / 0.7 ≈ 6.28 = T
#define PIL_EXAMPLE_IMAGE_TIME 4.4
#include "example_main.h"
#include <math.h>

typedef struct {
  pil_body *star;
  pil_body *planet;
} kp_ref;

static double sig_r(const pil_world *w, void *user) {
  const kp_ref *k = user;
  (void)w;

  return hypot(k->planet->pos.x - k->star->pos.x,
               k->planet->pos.y - k->star->pos.y);
}

static double sig_v(const pil_world *w, void *user) {
  const kp_ref *k = user;
  (void)w;

  return hypot(k->planet->vel.x, k->planet->vel.y);
}

static double sig_area(const pil_world *w, void *user) {
  const kp_ref *k = user;
  double x = k->planet->pos.x - k->star->pos.x;
  double y = k->planet->pos.y - k->star->pos.y;
  (void)w;

  return 0.5 * (x * k->planet->vel.y - y * k->planet->vel.x);
}

int main(int argc, char **argv) {
  const double GM = 1.0;
  const double e = 0.6;
  const double rp = 0.4;
  static kp_ref ref;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);

  double vp = sqrt(GM * (1.0 + e) / rp); // speed at perihelion
  double a = rp / (1.0 - e);
  double T = 2 * M_PI * sqrt(a * a * a / GM);

  pil_body *star =
      pil_add_body(w, &(pil_body_desc){.shape = pil_disc(0.07),
                                       .mass = 1.0,
                                       .pos = {0, 0},
                                       .flags = PIL_STATIC | PIL_HIDE_COM,
                                       .fill = 0.35});
  pil_body *planet = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(0.035),
                                                      .mass = 0.01,
                                                      .pos = {rp, 0},
                                                      .vel = {0, vp},
                                                      .fill = 0.9,
                                                      .flags = PIL_HIDE_COM});
  pil_add_gravity_pair(w, star, planet, GM / 1.0);

  pil_add_trace(w, pil_on(planet, 0, 0), 0.02, 600);
  pil_add_sweep_window(w, pil_on(star, 0, 0), pil_on(planet, 0, 0), 0.0, 0.5,
                       0.01);
  pil_add_sweep_window(w, pil_on(star, 0, 0), pil_on(planet, 0, 0),
                       0.5 * T - 0.25, 0.5 * T + 0.25, 0.01);
  pil_add_sweep(w, pil_on(star, 0, 0), pil_on(planet, 0, 0), 0.5, 0.01);
  pil_add_line(w, pil_on(star, 0, 0), pil_on(planet, 0, 0), 0);
  pil_add_vector(w, pil_on(planet, 0, 0), PIL_VEC_VELOCITY, 40.0, 0, "v");
  pil_add_vector(w, pil_on(planet, 0, 0), PIL_VEC_ACCELERATION, 8.0, 0, "a");
  pil_add_label(w, pil_at(0.0, -0.17), "S", 0, 0);
  pil_add_label(w, pil_at(rp + 0.06, 0.1), "perihelion", 2, 0);
  pil_add_label(w, pil_at(rp - 2 * a - 0.1, 0.1), "aphelion", -2, 0);

  ref.star = star;
  ref.planet = planet;
  pil_plot *pl =
      pil_add_plot(w, 0.66, 0.05, 0.99, 0.97, "t [s]", "r, v, dA/dt", T);
  pil_plot_add_signal(pl, "r", sig_r, &ref);
  pil_plot_add_signal(pl, "v", sig_v, &ref);
  pil_plot_add_signal(pl, "dA/dt", sig_area, &ref);
  pil_plots_record(w, 0.04, 5e-4);

  pil_view view;
  pil_view_fit_region(&view, -2.2, -1.0, 0.9, 1.0, 1500, 700, 0.0, 0.0, 0.64,
                      1.0, 12);
  int rc = example_main(w, &view, "Kepler orbit and equal areas", argc, argv,
                        "kepler.png");

  pil_world_free(w);

  return rc;
}
