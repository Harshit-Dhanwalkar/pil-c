/*
 * collision.c - oblique collision of two discs in lab frame (top) and in
 * centre-of-mass frame (bottom)
 * Momentum arrows m v on each disc; CM marker carries total momentum P = (m1 +
 * m2) v_cm: it moves in lab frame and is a point at rest, with P = 0, in CM
 * frame
 * Contact is frictionless with restitution e: only velocity component along
 * line of centres changes
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.5
#define PIL_EXAMPLE_IMAGE_TIME 2.1
#include "example_main.h"
#include <math.h>

// Callback for live label
static void theta_label_fn(const pil_world *w, void *user, char *buf, int n) {
  (void)w;
  pil_body **b = (pil_body **)user;
  double dx = b[1]->pos.x - b[0]->pos.x;
  double dy = b[1]->pos.y - b[0]->pos.y;
  double theta = atan2(dy, dx) * 180.0 / M_PI;
  if (theta < 0)
    theta += 360.0;
  snprintf(buf, n, "\xce\xb8 = %.0f\xc2\xb0", theta);
}

static pil_body *disc(pil_world *w, double r, double m, pil_vec2 p, pil_vec2 v,
                      double fill) {
  return pil_add_body(w, &(pil_body_desc){.shape = pil_disc(r),
                                          .mass = m,
                                          .pos = p,
                                          .vel = v,
                                          .fill = fill,
                                          .flags = PIL_HIDE_COM});
}

// Invisible walls for discs to bounce back into view instead of leaving
// viewport
static void add_walls(pil_world *w, pil_body *b) {
  // upper=0 means q >= limit, upper=1 means q <= limit
  pil_add_stop(w, b, PIL_COORD_X, -3.2, 0, 0.9); // Left wall
  pil_add_stop(w, b, PIL_COORD_X, 4.0, 1, 0.9);  // Right wall
  pil_add_stop(w, b, PIL_COORD_Y, -3.2, 0, 0.9); // Bottom wall
  pil_add_stop(w, b, PIL_COORD_Y, 3.6, 1, 0.9);  // Top wall
}

static void row(pil_world *w, double y, double vcm, double e) {
  const double r1 = 0.3;
  const double r2 = 0.45;
  const double m1 = 1.0;
  const double m2 = 2.0;
  const double u = 2.0;
  static pil_body *b[2][2];
  static int row_idx = 0;

  b[row_idx][0] =
      disc(w, r1, m1, (pil_vec2){-3.0, y}, (pil_vec2){u - vcm, 0}, 0.9);
  b[row_idx][1] =
      disc(w, r2, m2, (pil_vec2){0.5, y + 0.2}, (pil_vec2){-vcm, 0}, 0.6);

  pil_add_collision(w, b[row_idx][0], b[row_idx][1], r1, r2, e);

  // Dashed-line connecting centers
  pil_add_line(w, pil_on(b[row_idx][0], 0, 0), pil_on(b[row_idx][1], 0, 0), 1);

  // Angle of collision (theta) at center of first disc
  // It is measured from horizontal axis to line of centers
  pil_add_angle_between(w, pil_on(b[row_idx][0], 0, 0),
                        pil_on(b[row_idx][0], 1, 0),
                        pil_on(b[row_idx][1], 0, 0), 25, "\xce\xb8");

  for (int i = 0; i < 2; i++) {
    add_walls(w, b[row_idx][i]); // Keep them on screen
    pil_add_trace(w, pil_on(b[row_idx][i], 0, 0), 0.02, 2000);
    pil_add_vector(w, pil_on(b[row_idx][i], 0, 0), PIL_VEC_MOMENTUM, 35.0, 0,
                   i ? "p_2" : "p_1");
  }

  pil_add_live_label(w, pil_at(-3.4, y + 0.9), theta_label_fn, b[row_idx], 0,
                     0);

  row_idx++;
}

int main(int argc, char **argv) {
  const double e = 0.8;
  const double vcm = 2.0 * 1.0 / 3.0; // (m1 u) / (m1 + m2)
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);

  row(w, 1.9, 0.0, e);
  row(w, -1.9, vcm, e);

  pil_add_label(w, pil_at(-3.4, 3.4), "Lab frame", 2, 0);
  pil_add_label(w, pil_at(-3.4, 0.0), "CM frame (P = 0)", 2, 0);
  pil_add_label(w, pil_at(3.2, 3.4), "e = 0.8, m_2 = 2 m_1", 0, 0);

  pil_view view;
  pil_view_fit_region(&view, -3.6, -3.4, 4.2, 3.8, 1100, 1000, 0.0, 0.0, 1.0,
                      1.0, 12);
  int rc =
      example_main(w, &view, "Two-disc collision", argc, argv, "collision.png");

  pil_world_free(w);

  return rc;
}
