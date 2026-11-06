/*
 * impact.c - impact oscillator with a double pendulum hanging from it
 * Block m_1 slides on a rail between two stops (gap delta, restitution e) on a
 * spring k. A rod carries a hinge P from which disc m_2 (radius R, angle phi)
 * hangs by its rim; rod m_3 (length l, angle psi) is hinged at disc centre
 * Energy bar shows what each impact takes away (E_0, -dE, T, V)
 */

#include "example_main.h"
#include <math.h>

int main(int argc, char **argv) {
  const double wb = 0.5;   // block width
  const double hb = 0.4;   // block height
  const double stem = 0.8; // block's stem
  const double R = 0.35;
  const double l = 0.8;
  const double m1 = 1.0;
  const double m2 = 2.0;
  const double m3 = 1.0;
  const double k = 40.0;
  const double delta = 0.6;
  const double e = 0.6;
  const double wall = 1.5;
  const double x0 = -0.45;
  const double v0 = 1.8;
  const double phi0 = 0.35;
  const double psi0 = -0.5;
  const double yb = hb / 2; // block centre height

  pil_world *w = pil_world_new();

  pil_add_ground(w, (pil_vec2){1.2, hb}, (pil_vec2){-wall, hb}); // rail
  pil_add_ground(w, (pil_vec2){-wall, 0.8}, (pil_vec2){-wall, -0.05});

  // Stops: two small static posts, spring clears them above
  double xs = delta + wb / 2;
  for (int i = -1; i <= 1; i += 2) {
    pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.1, 0.2),
                                     .flags = PIL_STATIC | PIL_HIDE_COM,
                                     .pos = {i * (xs + 0.05), 0.1},
                                     .fill = 0.75});
  }

  // Block with its stem, one polygon (T upside down)
  pil_vec2 T[] = {{-wb / 2, 0}, {-0.025, 0}, {-0.025, -stem}, {0.025, -stem},
                  {0.025, 0},   {wb / 2, 0}, {wb / 2, hb},    {-wb / 2, hb}};
  pil_vec2 c1 = pil_poly_centroid(T, 8);
  pil_body *blk = pil_add_body(w, &(pil_body_desc){.shape = pil_poly(T, 8),
                                                   .mass = m1,
                                                   .pos = {x0, c1.y},
                                                   .vel = {v0, 0},
                                                   .flags = PIL_HIDE_COM,
                                                   .label = "m_1",
                                                   .label_dx = 0,
                                                   .label_dy = -6});
  pil_add_prismatic(w, pil_on_at(blk, x0, yb), pil_at(0, yb), (pil_vec2){1, 0});
  pil_add_stop(w, blk, PIL_COORD_X, delta, 1, e);
  pil_add_stop(w, blk, PIL_COORD_X, -delta, 0, e);
  pil_element *sp =
      pil_add_spring(w, pil_at(-wall, 0.32), pil_on_at(blk, x0 - wb / 2, 0.32),
                     k, wall - wb / 2);
  pil_element_set_label_at(sp, "k", 0, -18);

  // Disc hung from P by its rim
  pil_vec2 P = {x0, -stem};
  pil_vec2 C = {P.x + R * sin(phi0), P.y - R * cos(phi0)};
  pil_body *disc = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(R),
                                                    .mass = m2,
                                                    .pos = C,
                                                    .angle = phi0,
                                                    .vel = {v0, 0},
                                                    .label = "m_2",
                                                    .label_dx = 60,
                                                    .label_dy = -20});
  pil_add_revolute(w, pil_on(disc, 0, R), pil_on_at(blk, x0, -stem));

  // Rod m_3 hinged at disc centre
  pil_vec2 M = {C.x + l / 2 * sin(psi0), C.y - l / 2 * cos(psi0)};
  pil_body *rod = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(l, 0.05),
                                                   .mass = m3,
                                                   .pos = M,
                                                   .angle = psi0 - M_PI / 2,
                                                   .vel = {v0, 0},
                                                   .fill = 0.35,
                                                   .flags = PIL_HIDE_COM,
                                                   .label = "m_3",
                                                   .label_dx = -36,
                                                   .label_dy = -6});
  pil_add_revolute(w, pil_on(rod, -l / 2, 0), pil_on(disc, 0, 0));

  // Annotations
  pil_add_line(w, pil_at(0, hb), pil_at(0, hb + 0.35), 1);
  pil_add_dimension(w, pil_at(0, hb + 0.25), pil_on_at(blk, x0, hb + 0.25), 10,
                    "x");
  pil_add_dimension(w, pil_at(wb / 2, hb + 0.25),
                    pil_at(wb / 2 + delta, hb + 0.25), -10, "\xce\xb4");
  pil_add_line(w, pil_at(wb / 2, hb), pil_at(wb / 2, hb + 0.35), 1);
  pil_add_line(w, pil_at(wb / 2 + delta, hb), pil_at(wb / 2 + delta, hb + 0.35),
               1);
  pil_add_leader(w, pil_at(xs + 0.1, 0.12), 34, 26, "e");
  pil_add_angle_between(w, pil_on_at(blk, x0, -stem), pil_on_at(blk, x0, -1.5),
                        pil_on(disc, 0, 0), 38, "\xcf\x86");
  pil_add_angle_marker(w, pil_on(disc, 0, 0), -M_PI / 2, rod, 70, "\xcf\x88");
  pil_add_line(w, pil_on_at(blk, x0, -stem), pil_on(disc, 0, 0), 1);
  pil_add_label(w, pil_on_at(blk, x0, -stem), "P", -14, 6);
  pil_add_dimension(w, pil_on(disc, 0, 0), pil_on(disc, R, 0), 12, "R");
  pil_add_dimension(w, pil_on(rod, -l / 2, 0), pil_on(rod, l / 2, 0), -22, "l");
  pil_add_gravity_indicator(w, (pil_vec2){-1.3, -1.6});

  // Energy bar: V is measured from x = 0, \phi = \psi = 0
  double g = 9.81;
  double vref = g * (m1 * c1.y + m2 * (-stem - R) + m3 * (-stem - R - l / 2));
  pil_add_energy_bar(w, (pil_vec2){1.55, -1.7}, 26, 120, vref);

  pil_view view;
  pil_view_fit(&view, -1.75, -1.95, 2.25, 1.0, 900, 660, 24);
  int rc =
      example_main(w, &view, "Impact oscillator", argc, argv, "impact.png");

  pil_world_free(w);

  return rc;
}
