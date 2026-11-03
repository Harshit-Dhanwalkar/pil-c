/*
 * pulleys.c - block m_1 on a spring, tied over a fixed pulley to a hanging disc
 * m_2 (spring k_2 to ceiling), which carries a sprung mass m_3; two rod
 * pendulums swing from block
 * Everything is public-API
 */

#include "example_main.h"
#include <math.h>

int main(int argc, char **argv) {
  const double yc = 3.0;  // ceiling
  const double xw = 4.8;  // wall
  const double y1 = 2.2;  // block line
  const double R = 0.35;  // disc and
  const double rp = 0.18; // fixed-pulley radii
  pil_world *w = pil_world_new();

  // Hatching is drawn on right of a->b, so ceiling runs right to left to be
  // hatched above and wall runs bottom to top
  pil_add_ground(w, (pil_vec2){xw, yc}, (pil_vec2){0.1, yc});
  pil_add_ground(w, (pil_vec2){xw, 1.7}, (pil_vec2){xw, yc});

  // Fixed pulley: static disc; rope leaves its top horizontally and leaves its
  // left side vertically
  pil_vec2 Pc = {1.35 + rp, y1 - rp};
  pil_add_body(w, &(pil_body_desc){.shape = pil_disc(rp),
                                   .mass = 1,
                                   .pos = Pc,
                                   .flags = PIL_STATIC | PIL_HIDE_COM});

  // Block m_1 on a horizontal guide, spring k_1 to wall
  pil_body *m1 = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.5, 0.4),
                                                  .mass = 1.0,
                                                  .pos = {3.1, y1},
                                                  .label = "m_1",
                                                  .flags = PIL_HIDE_COM});
  pil_add_prismatic(w, pil_on(m1, 0, 0), pil_at(0, y1), (pil_vec2){1, 0});
  pil_element *k1 =
      pil_add_spring(w, pil_on(m1, 0.25, 0), pil_at(xw, y1), 25.0, 1.0);
  pil_element_set_label_at(k1, "k_1", 0, 22);

  // Hanging disc m_2, kept upright on a vertical guide
  pil_body *m2 = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(R),
                                                  .mass = 2.0,
                                                  .pos = {1.0, 1.1},
                                                  .label = "m_2",
                                                  .label_dx = 0,
                                                  .label_dy = 0});
  pil_add_prismatic(w, pil_on(m2, 0, 0), pil_at(1.0, 0), (pil_vec2){0, 1});
  pil_element *k2 =
      pil_add_spring(w, pil_at(0.9, yc), pil_on(m2, -0.1, R), 30.0, 1.1);
  pil_element_set_label_at(k2, "k_2", -26, 0);

  // Rope: block -> over pulley -> right edge of disc
  pil_vec2 c1 = {Pc.x, y1}, c2 = {Pc.x - rp, Pc.y};
  double len = (3.1 - 0.25 - c1.x) + (c2.y - 1.1);
  pil_add_pulley_rope(w, pil_on(m1, -0.25, 0), pil_on(m2, R, 0), c1, c2, len);

  // m_3 hangs from m_2 on spring k_3
  pil_body *m3 = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(0.15),
                                                  .mass = 0.6,
                                                  .pos = {1.0, 0.35},
                                                  .label = "m_3",
                                                  .label_dx = 26,
                                                  .label_dy = -6,
                                                  .flags = PIL_HIDE_COM});
  pil_add_prismatic(w, pil_on(m3, 0, 0), pil_at(1.0, 0), (pil_vec2){0, 1});
  pil_element *k3 =
      pil_add_spring(w, pil_on(m2, 0, -R), pil_on(m3, 0, 0.15), 60.0, 0.4);
  pil_element_set_label_at(k3, "k_3", 22, 0);

  // Two rod pendulums from block's base
  const double lA = 0.8;
  const double lB = 0.55;
  const double thA = 0.12;
  const double thB = -0.55;
  pil_vec2 pin = {3.1, y1 - 0.2};
  pil_body *pA = pil_add_body(
      w, &(pil_body_desc){.shape = pil_point_mass(0.06),
                          .mass = 0.4,
                          .pos = {pin.x + lA * sin(thA), pin.y - lA * cos(thA)},
                          .label = "m_A",
                          .label_dx = 20,
                          .label_dy = 6});
  pil_body *pB = pil_add_body(
      w, &(pil_body_desc){.shape = pil_point_mass(0.06),
                          .mass = 0.3,
                          .pos = {pin.x + lB * sin(thB), pin.y - lB * cos(thB)},
                          .label = "m_B",
                          .label_dx = -22,
                          .label_dy = -10});
  pil_add_distance(w, pil_on(m1, 0, -0.2), pil_on(pA, 0, 0), lA, 0);
  pil_add_distance(w, pil_on(m1, 0, -0.2), pil_on(pB, 0, 0), lB, 0);

  // Annotations
  pil_add_line(w, pil_on(m1, 0, -0.2), pil_at(pin.x, pin.y - lA - 0.12), 1);
  pil_add_angle_between(w, pil_on(m1, 0, -0.2), pil_at(pin.x, pin.y - 1),
                        pil_on(pA, 0, 0), 62, "\xce\xb8_A");
  pil_add_angle_between(w, pil_on(m1, 0, -0.2), pil_at(pin.x, pin.y - 1),
                        pil_on(pB, 0, 0), 36, "\xce\xb8_B");
  pil_add_dimension(w, pil_at(xw, y1 + 0.45), pil_on(m1, -0.25, 0.45), -26,
                    "x_1");
  pil_add_dimension(w, pil_at(0.55, yc), pil_on(m2, -0.45, 0), 0, "x_2");
  pil_add_dimension(w, pil_on(m2, -0.45, 0), pil_on(m3, -0.45, 0), 0, "x_3");
  pil_add_gravity_indicator(w, (pil_vec2){xw + 0.5, 2.8});

  pil_view view;
  pil_view_fit(&view, 0.0, -0.1, xw + 0.9, yc + 0.2, 900, 620, 24);
  int rc = example_main(w, &view, "Pulley system with pendulums", argc, argv,
                "pulleys.png");

  pil_world_free(w);

  return rc;
}
