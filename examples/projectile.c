/*
 * projectile.c - two shots from one point, with and without air drag
 * Same launch speed and angle
 * One ball flies in vacuum (parabola), other feels quadratic drag F = -c |v| v
 * and falls short and steeper
 * Paths are traced, velocity is drawn live with its x and y components and
 * dragged ball also shows its acceleration
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.6
#define PIL_EXAMPLE_IMAGE_TIME 1.5
#include "example_main.h"
#include <math.h>

typedef struct {
  pil_body *ball;
  double c; // drag coefficient, 1/m
} drag_ref;

static void drag_force(pil_world *w, double t, void *user) {
  const drag_ref *d = user;
  double sp =
      sqrt(d->ball->vel.x * d->ball->vel.x + d->ball->vel.y * d->ball->vel.y);
  (void)w;
  (void)t;

  pil_apply_force(d->ball, d->ball->pos,
                  (pil_vec2){-d->c * d->ball->mass * sp * d->ball->vel.x,
                             -d->c * d->ball->mass * sp * d->ball->vel.y});
}

int main(int argc, char **argv) {
  const double g = 9.81;
  const double R = 0.25;
  static drag_ref dr;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, -g);

  pil_add_ground(w, (pil_vec2){50, 0}, (pil_vec2){-2, 0});
  pil_vec2 v0 = {14.0, 16.0};

  pil_body *vac = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(R),
                                                   .mass = 1.0,
                                                   .pos = {0, R},
                                                   .vel = v0,
                                                   .fill = 0.9,
                                                   .flags = PIL_HIDE_COM});
  pil_body *air = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(R),
                                                   .mass = 1.0,
                                                   .pos = {0, R},
                                                   .vel = v0,
                                                   .fill = 0.55,
                                                   .flags = PIL_HIDE_COM});
  pil_add_floor_contact(w, vac, pil_at(0, 0), (pil_vec2){1, 0}, 1e5, 3e3, 0.8);
  pil_add_floor_contact(w, air, pil_at(0, 0), (pil_vec2){1, 0}, 1e5, 3e3, 0.8);
  dr.ball = air;
  dr.c = 0.02;
  pil_add_custom_force(w, drag_force, &dr);

  pil_add_trace(w, pil_on(vac, 0, 0), 0.02, 400);
  pil_add_trace(w, pil_on(air, 0, 0), 0.02, 400);
  pil_add_vector(w, pil_on(vac, 0, 0), PIL_VEC_VELOCITY, 9.0,
                 PIL_VEC_COMPONENTS, "v");
  pil_add_vector(w, pil_on(air, 0, 0), PIL_VEC_VELOCITY, 9.0, 0, "v");
  pil_add_vector(w, pil_on(air, 0, 0), PIL_VEC_ACCELERATION, 7.0, 0, "a");
  pil_add_label(w, pil_at(21, 14.2), "no drag: lighter ball", 2, 0);
  pil_add_label(w, pil_at(21, 8.2), "drag F = -c |v| v: darker ball", 2, 0);
  pil_add_gravity_indicator(w, (pil_vec2){40, 14});

  pil_view view;
  pil_view_fit_region(&view, -2, -1.2, 54, 14.5, 1800, 700, 0.0, 0.0, 1.0, 1.0,
                      14);
  int rc = example_main(w, &view, "Projectile with and without drag", argc,
                        argv, "projectile.png");

  pil_world_free(w);

  return rc;
}
