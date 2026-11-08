/*
 * loop.c - loop-the-loop
 * A small disc (m, radius r) is released at rest from height h on a ramp with
 * friction mu, runs through a rounded foot onto the floor and into a circular
 * loop of radius R (e = 0: no bounce). Where the normal force N reaches zero
 * the disc leaves the wall and flies freely until it meets it again
 * Public API only: floor contacts with ranges, circle contacts
 */

#include "example_main.h"
#include <math.h>

int main(int argc, char **argv) {
  const double alpha = 35.0 * M_PI / 180.0; /* ramp angle */
  const double h = 1.3;                     /* ramp height */
  const double Rf = 0.6;                    /* foot fillet */
  const double Rl = 0.55;                   /*  loop */
  const double r = 0.07;
  const double m = 1.0;
  const double mu = 0.5;
  const double k = 4e4;   /* contact stiffness  */
  const double c = 60.0;  /* damping, e = 0 */
  const double gap = 1.3; /* flat run between foot and loop */
  const double loop_span = 235.0 * M_PI / 180.0; /* loop arc, CCW from bottom */

  // Ramp tangent to fillet; ramp line passes through (0, h)
  double s = (h - Rf * (1 - cos(alpha))) / sin(alpha);
  double xf = s * cos(alpha) + Rf * sin(alpha);
  pil_vec2 T1 = {xf - Rf * sin(alpha), Rf - Rf * cos(alpha)};
  double xc = xf + gap;

  pil_world *w = pil_world_new();

  pil_add_ground(w, (pil_vec2){0, h}, T1);                       /* ramp */
  pil_add_ground(w, (pil_vec2){xf, 0}, (pil_vec2){xc + 1.0, 0}); /* floor */

  // Loop wall, a thin partial ring (215°, not a full circle)
  pil_body *ring = pil_add_body(
      w, &(pil_body_desc){.shape = pil_arc_band(Rl, Rl + 0.03, -M_PI / 2,
                                                -M_PI / 2 + loop_span),
                          .flags = PIL_STATIC | PIL_HIDE_COM,
                          .fill = 0.5});
  pil_vec2 o = pil_body_shape_origin(ring);
  pil_body_set_state(ring, (pil_vec2){xc - o.x, Rl - o.y}, 0, (pil_vec2){0, 0},
                     0);

  // Rounded foot joining ramp and floor
  pil_body *foot = pil_add_body(
      w, &(pil_body_desc){
             .shape = pil_arc_band(Rf, Rf + 0.03, -M_PI / 2 - alpha, -M_PI / 2),
             .flags = PIL_STATIC | PIL_HIDE_COM,
             .fill = 0.5});
  pil_vec2 of = pil_body_shape_origin(foot);
  pil_body_set_state(foot, (pil_vec2){xf - of.x, Rf - of.y}, 0,
                     (pil_vec2){0, 0}, 0);

  // Disc starts on ramp, one radius off surface
  pil_vec2 n = {sin(alpha), cos(alpha)};
  pil_vec2 start = {r * n.x, h + r * n.y};
  pil_body *disc =
      pil_add_body(w, &(pil_body_desc){.shape = pil_disc(r),
                                       .mass = m,
                                       .pos = start,
                                       .flags = PIL_SPOKE | PIL_HIDE_COM,
                                       .label = "m",
                                       .label_dx = 14,
                                       .label_dy = 16});

  pil_element *e;
  e = pil_add_floor_contact(w, disc, pil_at(T1.x, T1.y),
                            (pil_vec2){cos(alpha), -sin(alpha)}, k, c, mu);
  pil_floor_set_range(e, -s - 1.0, 0);
  e = pil_add_circle_contact(w, disc, pil_at(xf, Rf), Rf, 1, k, c, mu);
  pil_circle_contact_set_arc(e, -M_PI / 2 - alpha, alpha);
  e = pil_add_floor_contact(w, disc, pil_at(xf, 0), (pil_vec2){1, 0}, k, c, mu);
  pil_floor_set_range(e, 0, 1e9);
  // Loop wall starts at bottom and runs counter-clockwise  let open so disc can
  // enter along floor
  e = pil_add_circle_contact(w, disc, pil_at(xc, Rl), Rl, 1, k, c, mu);
  pil_circle_contact_set_arc(e, -M_PI / 2, loop_span);

  // Annotations
  pil_add_dimension(w, pil_at(0, 0), pil_at(0, h), 16, "h");
  pil_add_label(w,
                pil_at(T1.x - 0.5 * s * cos(alpha) * 0.5,
                       T1.y + 0.5 * s * sin(alpha) * 0.5),
                "\xce\xbc", 22, 12);
  pil_add_dimension(w, pil_at(xc, Rl),
                    pil_at(xc + Rl * cos(0.6), Rl + Rl * sin(0.6)), 0, "R");
  pil_add_label(w, pil_at(xc, Rl), "e = 0", -34, -30);
  pil_add_gravity_indicator(w, (pil_vec2){xc + 0.95, 1.1});

  pil_view view;
  pil_view_fit(&view, -0.4, -0.2, xc + 1.1, h + 0.35, 1000, 520, 24);
  int rc = example_main(w, &view, "Loop-the-loop", argc, argv, "loop.png");

  pil_world_free(w);

  return rc;
}
