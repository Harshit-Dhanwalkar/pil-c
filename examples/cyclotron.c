/*
 * cyclotron.c - charges in a uniform magnetic field (B out of page)
 * F = q v x B is always perpendicular to v, so speed is constant and path  is a
 * circle of radius r = m v / (q B), traversed with period T = 2 pi m / (q B)
 * that does not depend on speed: slow and fast positive charge return to start
 * together
 * Negative charge turns other way Velocity and force are drawn live on fast
 * positive charge
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 5.0
#include "example_main.h"
#include <math.h>

void pil_traces_clear(pil_world *w);

typedef struct {
  pil_world *w;
  double T;
  double last_t;
} trace_reset_user;

static void reset_traces_overlay(cairo_t *cr, const pil_view *view,
                                 void *user) {
  (void)cr;
  (void)view;
  trace_reset_user *u = (trace_reset_user *)user;

  double current_t = pil_world_time(u->w);

  int n_prev = (int)(u->last_t / u->T);
  int n_cur = (int)(current_t / u->T);
  if (n_cur > n_prev) {
    pil_traces_clear(u->w);
  }
  u->last_t = current_t;
}

int main(int argc, char **argv) {
  const double m = 1.0, q = 1.0, B = 1.0;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);

  pil_add_field_marks(w, (pil_vec2){-3.1, -3.5}, (pil_vec2){3.1, 3.0}, 0.55,
                      PIL_FIELD_B, 1, 0);

  struct {
    double v, q;
    double fill;
    const char *tag;
  } c[3] = {{0.8, q, 0.95, ""}, {1.6, q, 0.95, ""}, {1.2, -q, 0.2, ""}};
  pil_body *p[3];
  for (int i = 0; i < 3; i++) {
    p[i] = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.07),
                                            .mass = m,
                                            .pos = {0, 0},
                                            .vel = {c[i].v, 0},
                                            .fill = c[i].fill,
                                            .flags = PIL_HIDE_COM});
    pil_add_lorentz(w, p[i], c[i].q, 0, 0, B);
    pil_add_trace(w, pil_on(p[i], 0, 0), 0.02, 1000);
  }

  pil_add_vector(w, pil_on(p[1], 0, 0), PIL_VEC_VELOCITY, 40.0, 0, "v");
  pil_add_vector(w, pil_on(p[1], 0, 0), PIL_VEC_FORCE, 40.0, 0, "F");
  pil_add_label(w, pil_at(0.0, 0.25), "start", 0, 0);
  pil_add_label(w, pil_at(0.0, 3.3), "B out of the page,  T = 2 pi m / (q B)",
                0, 0);
  pil_add_label(w, pil_at(-2.9, 2.6), "r = m v / (q B)", 2, 0);

  // Setup trace reset every period
  static trace_reset_user tr_user;
  tr_user.w = w;
  tr_user.T = 2 * M_PI * m / (q * B); // Period = 2*pi
  tr_user.last_t = 0.0;
  pil_world_add_overlay(w, reset_traces_overlay, &tr_user);

  pil_view view;
  pil_view_fit_region(&view, -3.2, -3.6, 3.2, 3.6, 1000, 1000, 0.0, 0.0, 1.0,
                      1.0, 12);
  int rc = example_main(w, &view, "Charges in a magnetic field", argc, argv,
                        "cyclotron.png");

  pil_world_free(w);

  return rc;
}
