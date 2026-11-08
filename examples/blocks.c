/*
 * blocks.c - two blocks on springs, each carrying a rod pendulum
 * Wall -- k -- m_A -- k -- m_B on a frictionless floor; pendulums l_1, l_2 hang
 * through floor line
 * Everything is public-API
 */

#include "example_main.h"
#include <math.h>

int main(int argc, char **argv) {
  const double hb = 0.2;  // block half-height
  const double wb = 0.55; // block width
  const double l1 = 0.75;
  const double l2 = 0.95;
  // block centres (initial)
  const double xa0 = 1.0;
  const double xb0 = 2.75;

  pil_world *w = pil_world_new();

  pil_add_ground(w, (pil_vec2){-0.5, 0}, (pil_vec2){4.2, 0});
  pil_add_ground(w, (pil_vec2){0, 0.8}, (pil_vec2){0, -0.05}); // wall

  pil_body *A = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(wb, 2 * hb),
                                                 .mass = 1.5,
                                                 .pos = {xa0, hb},
                                                 .label = "m_A",
                                                 .label_dy = 0});
  pil_body *B = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(wb, 2 * hb),
                                                 .mass = 1.0,
                                                 .pos = {xb0, hb},
                                                 .label = "m_B"});
  pil_add_prismatic(w, pil_on(A, 0, 0), pil_at(0, hb), (pil_vec2){1, 0});
  pil_add_prismatic(w, pil_on(B, 0, 0), pil_at(0, hb), (pil_vec2){1, 0});

  pil_element *s0 =
      pil_add_spring(w, pil_at(0, hb), pil_on(A, -wb / 2, 0), 12.0, 0.6);
  pil_element *s1 =
      pil_add_spring(w, pil_on(A, wb / 2, 0), pil_on(B, -wb / 2, 0), 12.0, 0.8);
  pil_element_set_label_at(s0, "k", 0, 22);
  pil_element_set_label_at(s1, "k", 0, 22);

  // Pendulums: rigid massless rods from a pin at block's base
  double th1 = 0.35, th2 = -0.25;
  pil_body *p1 = pil_add_body(
      w, &(pil_body_desc){.shape = pil_point_mass(0.06),
                          .mass = 0.5,
                          .pos = {xa0 + l1 * sin(th1), -l1 * cos(th1)},
                          .label = "m_1",
                          .label_dx = 18,
                          .label_dy = -14});
  pil_body *p2 = pil_add_body(
      w, &(pil_body_desc){.shape = pil_point_mass(0.06),
                          .mass = 0.4,
                          .pos = {xb0 + l2 * sin(th2), -l2 * cos(th2)},
                          .label = "m_2",
                          .label_dx = 18,
                          .label_dy = -14});
  pil_add_distance(w, pil_on(A, 0, -hb), pil_on(p1, 0, 0), l1, 0);
  pil_add_distance(w, pil_on(B, 0, -hb), pil_on(p2, 0, 0), l2, 0);

  // Annotations
  pil_add_dimension(w, pil_at(0, 2 * hb), pil_on(A, -wb / 2, hb), 34, "x_A");
  pil_add_dimension(w, pil_at(0, 2 * hb), pil_on(B, -wb / 2, hb), 64, "x_B");
  pil_add_line(w, pil_on(A, 0, -hb), pil_at(xa0, -l1 - 0.08), 1);
  pil_add_line(w, pil_on(B, 0, -hb), pil_at(xb0, -l2 - 0.08), 1);
  pil_add_angle_between(w, pil_on(A, 0, -hb), pil_at(xa0, -1), pil_on(p1, 0, 0),
                        36, "\xce\xb8_1");
  pil_add_angle_between(w, pil_on(B, 0, -hb), pil_at(xb0, -1), pil_on(p2, 0, 0),
                        36, "\xce\xb8_2");
  pil_add_label(w, pil_on(p1, 0, 0), "l_1", 14, 44);
  pil_add_label(w, pil_on(p2, 0, 0), "l_2", 14, 52);

  pil_view view;
  pil_view_fit(&view, -0.4, -1.25, 4.1, 1.35, 900, 480, 24);
  int rc = example_main(w, &view, "Two blocks with pendulums", argc, argv,
                        "blocks.png");

  pil_world_free(w);

  return rc;
}
