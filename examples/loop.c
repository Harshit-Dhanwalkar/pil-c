/*
 * loop.c - loop-the-loop
 * A small disc (m, radius r) starts on a HORIZONTAL PLANK at height h with a
 * small initial push, rolls off end, down a slope with friction mu, through a
 * rounded foot onto floor and into a partial circular loop of radius R (e = 0:
 * no bounce). Where normal force N reaches zero disc leaves wall and flies
 * freely until it meets it again
 * Public API only: floor contacts with ranges, circle contacts
 */

#include "example_main.h"
#include <math.h>

int main(int argc, char **argv) {
  const double alpha = 35.0 * M_PI / 180.0; /* slope angle (after plank) */
  const double h = 1.3;                     /* plank / slope height */
  const double Rf = 0.6;                    /* foot fillet */
  const double Rl = 0.55;                   /* loop radius */
  const double r = 0.07;
  const double m = 1.0;
  const double mu = 0.5;
  const double k = 4e4;   /* contact stiffness  */
  const double c = 60.0;  /* damping, e = 0 */
  const double gap = 1.3; /* flat run between foot and loop */
  const double loop_span = 270.0 * M_PI / 180.0; /* loop arc, CCW from bottom */
  const double plank_len = 0.8;                  /* horizontal plank at top */
  const double v0 = 3.0; /* initial push so disc traverses plank */

  // Slope starts at end of plank, (plank_len, h), and runs down at
  // angle alpha to top of fillet
  // Slope is tangent to fillet
  double s = (h - Rf * (1 - cos(alpha))) / sin(alpha);
  double xf = plank_len + s * cos(alpha) + Rf * sin(alpha);
  pil_vec2 T1 = {xf - Rf * sin(alpha), Rf - Rf * cos(alpha)};
  pil_vec2 T_plank = {plank_len, h}; /* end of plank = top of slope */
  double xc = xf + gap;

  pil_world *w = pil_world_new();

  pil_add_ground(w, (pil_vec2){0, h}, T_plank);                  /* plank */
  pil_add_ground(w, T_plank, T1);                                /* slope */
  pil_add_ground(w, (pil_vec2){xf, 0}, (pil_vec2){xc + 1.0, 0}); /* floor */

  // Loop wall, a thin partial ring
  pil_body *ring = pil_add_body(
      w, &(pil_body_desc){.shape = pil_arc_band(Rl, Rl + 0.03, -M_PI / 2,
                                                -M_PI / 2 + loop_span),
                          .flags = PIL_STATIC | PIL_HIDE_COM,
                          .fill = 0.5});
  pil_vec2 o = pil_body_shape_origin(ring);
  pil_body_set_state(ring, (pil_vec2){xc - o.x, Rl - o.y}, 0, (pil_vec2){0, 0},
                     0);

  // Rounded foot joining slope and floor
  pil_body *foot = pil_add_body(
      w, &(pil_body_desc){
             .shape = pil_arc_band(Rf, Rf + 0.03, -M_PI / 2 - alpha, -M_PI / 2),
             .flags = PIL_STATIC | PIL_HIDE_COM,
             .fill = 0.5});
  pil_vec2 of = pil_body_shape_origin(foot);
  pil_body_set_state(foot, (pil_vec2){xf - of.x, Rf - of.y}, 0,
                     (pil_vec2){0, 0}, 0);

  // Disc starts on horizontal plank, one radius above it, with initial push
  pil_vec2 start = {0.15, h + r};
  pil_body *disc =
      pil_add_body(w, &(pil_body_desc){.shape = pil_disc(r),
                                       .mass = m,
                                       .pos = start,
                                       .flags = PIL_SPOKE | PIL_HIDE_COM,
                                       .label = "m",
                                       .label_dx = 14,
                                       .label_dy = 16});
  pil_body_set_state(disc, start, 0, (pil_vec2){v0, 0}, 0);

  pil_element *e;
  // Horizontal plank contact (disc on top, moving right)
  e = pil_add_floor_contact(w, disc, pil_at(0, h), (pil_vec2){1, 0}, k, c, mu);
  pil_floor_set_range(e, -0.1, plank_len + 0.05);
  // Slope contact (from T1 up-left to plank end)
  e = pil_add_floor_contact(w, disc, pil_at(T1.x, T1.y),
                            (pil_vec2){cos(alpha), -sin(alpha)}, k, c, mu);
  pil_floor_set_range(e, -s - 1.0, 0);
  // Fillet
  e = pil_add_circle_contact(w, disc, pil_at(xf, Rf), Rf, 1, k, c, mu);
  pil_circle_contact_set_arc(e, -M_PI / 2 - alpha, alpha);
  // Floor
  e = pil_add_floor_contact(w, disc, pil_at(xf, 0), (pil_vec2){1, 0}, k, c, mu);
  pil_floor_set_range(e, 0, 1e9);
  // Loop wall (partial arc, starts at bottom and runs CCW)
  e = pil_add_circle_contact(w, disc, pil_at(xc, Rl), Rl, 1, k, c, mu);
  pil_circle_contact_set_arc(e, -M_PI / 2, loop_span);

  // Annotations
  pil_add_dimension(w, pil_at(0, 0), pil_at(0, h), 16, "h");
  pil_add_label(
      w, pil_at(T1.x - 0.25 * s * cos(alpha), T1.y + 0.25 * s * sin(alpha)),
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
