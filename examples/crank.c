/*
 * crank.c - a torsion-sprung crank disc drives two rows of sliding masses
 * Disc M (spring kappa) -> crank pin P -> rod l -> block m in top channel;
 * springs chain blocks, and diagonal springs couple top row to  bottom row
 * Everything is public-API
 */

#include "example_main.h"
#include <math.h>

// Block sliding on a horizontal line y = cy, centred at x
static pil_body *slider(pil_world *w, double x, double cy, const char *label) {
  pil_body *b = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.5, 0.5),
                                                 .mass = 1.0,
                                                 .pos = {x, cy},
                                                 .label = label,
                                                 .flags = PIL_HIDE_COM});
  pil_add_prismatic(w, pil_on(b, 0, 0), pil_at(0, cy), (pil_vec2){1, 0});

  return b;
}

int main(int argc, char **argv) {
  const double xr = 6.0;
  const double ytop = 1.95;
  const double ybot = 0.5;
  pil_world *w = pil_world_new();

  // Channels: hatched walls above and below each row
  pil_add_ground(w, (pil_vec2){1.7, ytop + 0.32}, (pil_vec2){xr, ytop + 0.32});
  pil_add_ground(w, (pil_vec2){xr, ytop - 0.32}, (pil_vec2){1.7, ytop - 0.32});
  pil_add_ground(w, (pil_vec2){0.4, ybot + 0.32}, (pil_vec2){xr, ybot + 0.32});
  pil_add_ground(w, (pil_vec2){xr, ybot - 0.32}, (pil_vec2){0.4, ybot - 0.32});

  pil_add_ground(w, (pil_vec2){xr, ybot - 0.32}, (pil_vec2){xr, ytop + 0.32});

  // disc with crank pin at radius d, torsion spring kappa to world
  const double d = 0.4;
  const double Rd = 0.62;
  pil_vec2 O = {0.85, ytop};
  pil_body *disc = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(Rd),
                                                    .mass = 1.5,
                                                    .pos = O,
                                                    .flags = PIL_HIDE_COM,
                                                    .label = "M",
                                                    .label_dx = -42,
                                                    .label_dy = -34});
  pil_add_revolute(w, pil_on(disc, 0, 0), pil_at(O.x, O.y));
  pil_add_torsion_spring(w, disc, NULL, 3.0, 0.9, 0.04);
  pil_vec2 P = {O.x + d * cos(0.5), O.y + d * sin(0.5)};

  pil_body *A = slider(w, 2.5, ytop, "m");
  pil_body *B = slider(w, 4.3, ytop, "m");
  pil_body *C = slider(w, 2.9, ybot, "m");
  pil_body *D = slider(w, 4.9, ybot, "m");

  // Connecting rod: light bar from crank pin to left face of A
  pil_vec2 Q = {2.25, ytop};
  double L = hypot(Q.x - P.x, Q.y - P.y);
  pil_body *rod = pil_add_body(
      w, &(pil_body_desc){.shape = pil_rect(L, 0.04),
                          .mass = 0.2,
                          .pos = {(P.x + Q.x) / 2, (P.y + Q.y) / 2},
                          .angle = atan2(Q.y - P.y, Q.x - P.x),
                          .fill = 0.1,
                          .flags = PIL_HIDE_COM});
  pil_add_revolute(w, pil_on_at(disc, P.x, P.y), pil_on_at(rod, P.x, P.y));
  pil_add_revolute(w, pil_on_at(rod, Q.x, Q.y), pil_on_at(A, Q.x, Q.y));

  // Top row
  pil_element *s;
  s = pil_add_spring(w, pil_on(A, 0.25, 0), pil_on(B, -0.25, 0), 20.0, 1.0);
  pil_element_set_label_at(s, "k_3", 0, 20);
  s = pil_add_spring(w, pil_on(B, 0.25, 0), pil_at(xr, ytop), 25.0, 1.3);
  pil_element_set_label_at(s, "k_2", 0, 20);
  // Bottom row
  s = pil_add_spring(w, pil_at(0.4, ybot), pil_on(C, -0.25, 0), 22.0, 1.8);
  pil_element_set_label_at(s, "k_1", 0, 20);
  s = pil_add_spring(w, pil_on(C, 0.25, 0), pil_on(D, -0.25, 0), 22.0, 1.4);
  pil_element_set_label_at(s, "k_1", 0, 20);
  s = pil_add_spring(w, pil_on(D, 0.25, 0), pil_at(xr, ybot), 22.0, 0.7);
  pil_element_set_label_at(s, "k_1", 0, 20);
  // Diagonal couplings
  s = pil_add_spring(w, pil_on(A, 0.1, -0.25), pil_on(C, 0.0, 0.25), 18.0,
                     PIL_REST_INITIAL);
  pil_element_set_label_at(s, "k_2", 22, 6);
  s = pil_add_spring(w, pil_on(B, 0.1, -0.25), pil_on(D, 0.0, 0.25), 18.0,
                     PIL_REST_INITIAL);
  pil_element_set_label_at(s, "k_2", 22, 6);

  // Annotations
  pil_add_dimension(w, pil_at(xr, ytop + 0.32), pil_on(A, -0.25, 0.32), -22,
                    "x_1");
  pil_add_dimension(w, pil_at(xr, ytop + 0.32), pil_on(B, -0.25, 0.32), -48,
                    "x_2");
  pil_add_dimension(w, pil_at(xr, ybot - 0.32), pil_on(D, -0.25, -0.32), 22,
                    "x_3");
  pil_add_dimension(w, pil_at(xr, ybot - 0.32), pil_on(C, -0.25, -0.32), 48,
                    "x_4");
  pil_add_dimension(w, pil_at(O.x, O.y), pil_on_at(disc, P.x, P.y), 20, "d");
  pil_add_dimension(w, pil_on_at(rod, P.x, P.y), pil_on_at(rod, Q.x, Q.y), -26,
                    "\xe2\x84\x93");
  pil_add_label(w, pil_at(O.x, O.y), "\xce\xba", 0, 0);
  pil_add_gravity_indicator(w, (pil_vec2){xr + 0.55, 2.7});

  pil_view view;
  pil_view_fit(&view, -0.1, -0.35, xr + 0.9, 3.0, 1000, 620, 24);
  int rc = example_main(w, &view, "Crank-driven coupled sliders", argc, argv,
                        "crank.png");

  pil_world_free(w);

  return rc;
}
