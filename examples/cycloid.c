/*
 * cycloid.c - rolling wheel and paths of points fixed to it
 * A disc rolls without slipping on ground at constant speed
 * Points on rim, inside it and outside it (a flange) trace a cycloid, a curtate
 * and a prolate cycloid Rim point's velocity is drawn live with its components:
 * it vanishes whenever that point touches ground
 * The rolling-without-slipping constraint links angle to x, so reversing
 * x also reverses theta
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.7
#define PIL_EXAMPLE_IMAGE_TIME 3.2
#include "example_main.h"
#include <math.h>

static const double R = 0.5;
static const double PEAK_SPEED = 2.0;
static const double TRACK_LEFT = 0.0;
static const double TRACK_RIGHT = 11.5;

// x(t) = TRACK_LEFT + A*(1 - cos(w t)),  A*(w) = PEAK_SPEED
// Oscillates between TRACK_LEFT and TRACK_RIGHT with peak speed PEAK_SPEED
static double wheel_x_drive(double t, void *user) {
  (void)user;
  double A = 0.5 * (TRACK_RIGHT - TRACK_LEFT);
  double w = PEAK_SPEED / A;
  return TRACK_LEFT + A * (1.0 - cos(w * t));
}

int main(int argc, char **argv) {
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);

  pil_add_ground(w, (pil_vec2){12, 0}, (pil_vec2){-1.5, 0});

  pil_body *wheel = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(R),
                                                     .mass = 1.0,
                                                     .pos = {TRACK_LEFT, R},
                                                     .vel = {0, 0},
                                                     .omega = 0,
                                                     .flags = PIL_SPOKE});

  // No-slip on the ground line:  x + R*theta = const  (phase = 0 here)
  pil_add_rolling(w, wheel, (pil_vec2){0, 0}, R, (pil_vec2){0, 0},
                  (pil_vec2){1, 0});

  // Drive x
  pil_add_coord_drive(w, wheel, PIL_COORD_X, wheel_x_drive, NULL);

  // Tracer curves
  pil_add_trace(w, pil_on(wheel, R, 0), 0.01, 4000);
  pil_add_trace(w, pil_on(wheel, 0.5 * R, 0), 0.01, 4000);
  pil_add_trace(w, pil_on(wheel, 1.5 * R, 0), 0.01, 4000);

  pil_add_line(w, pil_on(wheel, 0, 0), pil_on(wheel, 1.5 * R, 0), 1);
  pil_add_vector(w, pil_on(wheel, R, 0), PIL_VEC_VELOCITY, 40.0,
                 PIL_VEC_COMPONENTS, "v");
  pil_add_arrow(w, pil_on(wheel, 0, 0), 60, 0, "v_C");
  pil_add_label(w, pil_at(-1.3, 1.5), "rim: cycloid", 2, 0);
  pil_add_label(w, pil_at(-1.3, 1.25), "inside: curtate", 2, 0);
  pil_add_label(w, pil_at(-1.3, 1.0), "outside: prolate", 2, 0);

  pil_view view;
  pil_view_fit_region(&view, -1.5, -0.4, 12, 2.0, 1500, 500, 0.0, 0.0, 1.0, 1.0,
                      14);
  int rc = example_main(w, &view, "Rolling wheel and its cycloids", argc, argv,
                        "cycloid.png");

  pil_world_free(w);

  return rc;
}
