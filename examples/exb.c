/*
 * exb.c - E x B drift of a charge in crossed uniform fields
 * E points up, B out of page
 * Charge is released at rest: it is accelerated along E, curved by B and motion
 * is a cycloid Averaged over a turn charge drifts sideways at v_d = E / B,
 * perpendicular to both fields (same for any charge or mass)
 * Plot shows v_x oscillating about E / B
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 9.0
#include "example_main.h"
#include <math.h>

static double sig_vx(const pil_world *w, void *user) {
  (void)w;
  return ((pil_body *)user)->vel.x;
}

static double sig_vd(const pil_world *w, void *user) {
  (void)w;
  (void)user;
  return 0.5; // E / B
}

static void loop_world_overlay(cairo_t *cr, const pil_view *view, void *user) {
  (void)cr;
  (void)view;
  pil_world *w = (pil_world *)user;
  if (pil_world_time(w) >= 12.0) {
    pil_world_reset(w);
  }
}

int main(int argc, char **argv) {
  const double E = 0.5, B = 1.0;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);

  pil_add_field_marks(w, (pil_vec2){-0.5, -0.4}, (pil_vec2){6.5, 1.5}, 0.5,
                      PIL_FIELD_B, 1, 0);
  pil_add_field_marks(w, (pil_vec2){-0.5, -0.4}, (pil_vec2){6.5, 1.5}, 1.25,
                      PIL_FIELD_E, 0, 1);

  pil_body *p = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.06),
                                                 .mass = 1.0,
                                                 .pos = {0, 0.1},
                                                 .fill = 0.95,
                                                 .flags = PIL_HIDE_COM});
  pil_add_lorentz(w, p, 1.0, 0, E, B);
  pil_add_trace(w, pil_on(p, 0, 0), 0.02, 1500);
  pil_add_vector(w, pil_on(p, 0, 0), PIL_VEC_VELOCITY, 50.0, PIL_VEC_COMPONENTS,
                 "v");
  pil_add_label(w, pil_at(0.0, -0.25), "start at rest", 2, 0);
  pil_add_label(w, pil_at(5.5, 1.7), "drift v_d = E / B", 0, 0);

  pil_plot *pl = pil_add_plot(w, 0.08, 0.05, 0.98, 0.38, "t [s]", "v_x", 12.0);
  pil_plot_add_signal(pl, "v_x", sig_vx, p);
  pil_plot_add_signal(pl, "E/B", sig_vd, NULL);
  pil_plots_record(w, 0.05, 5e-4);

  // Register loop overlay to reset world when it reaches plot's span
  pil_world_add_overlay(w, loop_world_overlay, w);

  pil_view view;
  pil_view_fit_region(&view, -0.7, -0.5, 11.2, 1.9, 1500, 800, 0.0, 0.4, 1.0,
                      1.0, 12);
  int rc = example_main(w, &view, "E x B drift", argc, argv, "exb.png");

  pil_world_free(w);

  return rc;
}
