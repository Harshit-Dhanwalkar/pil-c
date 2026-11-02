/*
 * frame.c - a sliding frame M holding two sprung blocks, with a two-link
 * pendulum (torsion knee kappa) on its left and a sprung link on its right,
 * pulled toward a wall by k_4. Everything is public-API; every point is given
 * in world coordinates and pil_on_at() turns it into a body anchor
 */

#include "example_main.h"
#include <math.h>

/* thick black link between world points p and q, mass m */
static pil_body *link(pil_world *w, pil_vec2 p, pil_vec2 q, double m,
                      const char *label, double ldx, double ldy) {
  double len = hypot(q.x - p.x, q.y - p.y);

  return pil_add_body(
      w, &(pil_body_desc){.shape = pil_rect(len, 0.055),
                          .mass = m,
                          .pos = {(p.x + q.x) / 2, (p.y + q.y) / 2},
                          .angle = atan2(q.y - p.y, q.x - p.x),
                          .fill = 0.08,
                          .flags = PIL_HIDE_COM,
                          .label = label,
                          .label_dx = ldx,
                          .label_dy = ldy});
}

int main(int argc, char **argv) {
  const double wall_x = 6.0;
  pil_world *w = pil_world_new();

  pil_add_ground(w, (pil_vec2){0.1, 0}, (pil_vec2){wall_x, 0});
  pil_add_ground(w, (pil_vec2){wall_x, 0}, (pil_vec2){wall_x, 2.7});

  // U-shaped frame, authored by its outline
  pil_vec2 U[8] = {{2.0, 0},  {2.14, 0}, {2.14, 2.06}, {2.86, 2.06},
                   {2.86, 0}, {3.0, 0},  {3.0, 2.2},   {2.0, 2.2}};
  pil_body *M = pil_add_body(w, &(pil_body_desc){.shape = pil_poly(U, 8),
                                                 .mass = 6.0,
                                                 .pos = pil_poly_centroid(U, 8),
                                                 .flags = PIL_HIDE_COM});
  pil_add_prismatic(w, pil_on(M, 0, 0), pil_at(0, M->pos.y), (pil_vec2){1, 0});

  // 2 blocks guided vertically inside frame
  pil_body *m1 = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.70, 0.5),
                                                  .mass = 1.0,
                                                  .pos = {2.5, 1.35},
                                                  .label = "m"});
  pil_body *m2 = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.70, 0.5),
                                                  .mass = 1.0,
                                                  .pos = {2.5, 0.65},
                                                  .label = "m"});
  pil_add_prismatic(w, pil_on(m1, 0, 0), pil_on_at(M, 2.5, 1.0),
                    (pil_vec2){0, 1});
  pil_add_prismatic(w, pil_on(m2, 0, 0), pil_on_at(M, 2.5, 1.0),
                    (pil_vec2){0, 1});
  pil_element *k1a = pil_add_spring(w, pil_on_at(M, 2.5, 2.06),
                                    pil_on(m1, 0, 0.25), 45.0, 0.62);
  pil_element *k1b =
      pil_add_spring(w, pil_on(m1, 0, -0.25), pil_on(m2, 0, 0.25), 80.0, 0.25);
  pil_element_set_label_at(k1a, "k_1", 22, 0);
  pil_element_set_label_at(k1b, "k_1", 22, 0);

  // Left: two links, hinge on frame, knee with torsion spring kappa
  pil_vec2 H1 = {2.0, 1.45}, K = {1.2, 0.85}, T = {0.45, 0.2};
  pil_body *r1 = link(w, H1, K, 0.6, "m", -14, 24);
  pil_body *r2 = link(w, K, T, 0.6, "m", -22, 22);
  pil_add_revolute(w, pil_on_at(M, H1.x, H1.y), pil_on_at(r1, H1.x, H1.y));
  pil_add_revolute(w, pil_on_at(r1, K.x, K.y), pil_on_at(r2, K.x, K.y));
  pil_add_torsion_spring(w, r1, r2, 2.0, r1->angle - r2->angle, 0.02);
  pil_element *k2 = pil_add_spring(w, pil_on_at(r1, K.x, K.y),
                                   pil_on_at(M, 2.0, 0.5), 30.0, 0.45);
  pil_element_set_label_at(k2, "k_2", -10, -26);
  pil_add_label(w, pil_on_at(r1, K.x, K.y), "\xce\xba", -26, 26);

  // Right: one link hinged on frame, tip sprung to frame foot
  pil_vec2 H2 = {3.0, 1.45}, T2 = {3.3, 0.85};
  pil_body *r3 = link(w, H2, T2, 0.6, "m", 24, 12);
  pil_add_revolute(w, pil_on_at(M, H2.x, H2.y), pil_on_at(r3, H2.x, H2.y));
  pil_element *k3 = pil_add_spring(w, pil_on_at(r3, T2.x, T2.y),
                                   pil_on_at(M, 3.0, 0.45), 30.0, 0.4);
  pil_element_set_label_at(k3, "k_3", 14, -24);

  // k_4 from right hinge to wall
  pil_element *k4 = pil_add_spring(w, pil_on_at(M, H2.x, H2.y),
                                   pil_at(wall_x, H2.y), 12.0, 2.6);
  pil_element_set_label_at(k4, "k_4", 0, 26);

  pil_add_leader(w, pil_on_at(M, 2.75, 2.2), 26, 26, "M");
  pil_add_gravity_indicator(w, (pil_vec2){wall_x + 0.55, 2.6});

  pil_view view;
  pil_view_fit(&view, 0.0, -0.1, wall_x + 0.8, 2.85, 1000, 590, 24);
  int rc = example_main(w, &view, "Sliding frame with sprung links", argc, argv,
                        "frame.png");

  pil_world_free(w);

  return rc;
}
