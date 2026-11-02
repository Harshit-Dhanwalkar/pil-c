/*
 * ring.c - plan view: a hub M with a torsion spring, ringed by point masses
 * joined by springs (g points into page, so nothing falls)
 * Everything is public-API
 */

#include "example_main.h"
#include <math.h>

#define N 12

int main(int argc, char **argv) {
  const double Rh = 0.45;
  const double Rin = 0.9;
  const double Rout = 1.35;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);

  // Hub: regular 12-gon, pinned at origin, torsion spring to world
  pil_vec2 hub[N];
  for (int i = 0; i < N; i++) {
    hub[i] = (pil_vec2){Rh * cos(2 * M_PI * i / N), Rh * sin(2 * M_PI * i / N)};
  }
  pil_body *M = pil_add_body(w, &(pil_body_desc){.shape = pil_poly(hub, N),
                                                 .mass = 2.0,
                                                 .label = "M",
                                                 .label_dx = 0,
                                                 .label_dy = 16,
                                                 .flags = PIL_HIDE_COM});
  pil_add_revolute(w, pil_on(M, 0, 0), pil_at(0, 0));
  pil_add_torsion_spring(w, M, NULL, 4.0, 0.7, 0.05);
  // hub starts wound up by 0.7 rad: rest angle is offset from 0

  // masses alternate between an inner and an outer ring
  pil_body *m[N];
  for (int i = 0; i < N; i++) {
    double a = 2 * M_PI * (i + 0.5) / N, R = (i % 2) ? Rout : Rin;
    m[i] = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.055),
                                            .mass = 0.3,
                                            .pos = {R * cos(a), R * sin(a)},
                                            .flags = PIL_HIDE_COM});
  }

  for (int i = 0; i < N; i += 2) {
    // inner masses ride on rigid links from hub
    // springs tie them to each other and to free outer masses between them
    pil_add_distance(w, pil_on(M, hub[i].x, hub[i].y), pil_on(m[i], 0, 0),
                     PIL_LENGTH_INITIAL, 0);

    pil_element *sp1 =
        pil_add_spring(w, pil_on(m[i], 0, 0), pil_on(m[(i + 2) % N], 0, 0),
                       25.0, PIL_REST_INITIAL);
    pil_element_set_label(sp1, "k");

    pil_element *sp2 = pil_add_spring(
        w, pil_on(m[i], 0, 0), pil_on(m[i + 1], 0, 0), 25.0, PIL_REST_INITIAL);
    pil_element_set_label(sp2, "k");

    pil_element *sp3 =
        pil_add_spring(w, pil_on(m[i + 1], 0, 0), pil_on(m[(i + 2) % N], 0, 0),
                       25.0, PIL_REST_INITIAL);
    pil_element_set_label(sp3, "k");
  }

  pil_add_label(w, pil_on(m[0], 0, 0), "m", 12, 10);
  // pil_add_label(w, pil_at(0, 0), "\xce\xba", 0, -16);
  pil_add_gravity_into_page(w, (pil_vec2){1.75, 1.45});

  pil_view view;
  pil_view_fit(&view, -2.0, -1.7, 2.0, 1.7, 900, 700, 24);
  int rc = example_main(w, &view, "Ring of masses around a hub", argc, argv,
                        "ring.png");

  pil_world_free(w);

  return rc;
}
