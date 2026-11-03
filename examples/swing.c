/*
 * swing.c - a swinging Atwood machine that drives a shuttle
 *
 * A ball m hangs on a rope over two pulleys P1, P2 and down to a block M that
 * rests on a spring k
 * A belt from P2 turns two pinch rollers whose grip moves a capsule m_c back
 * and forth
 * Pulleys, rollers and capsule are real bodies tied to rope speed by
 * pil_add_link(), so their inertia is part of dynamics; spring is one that
 * makes bounce exactly twice as fast as swing
 * Panel below figure plots energy
 * Everything is public-API
 */

#include "example_main.h"
#include <math.h>

typedef struct {
  pil_body *m;
  pil_vec2 pivot;
  double mass, g, E0;
} swing_ref;

// Energy of swing alone: 1/2 m (r \theta')^2 + m g r (1 - \cos \theta)
static double e_swing(const pil_world *w, void *user) {
  const swing_ref *s = user;
  pil_vec2 d = {s->m->pos.x - s->pivot.x, s->m->pos.y - s->pivot.y};
  double r = hypot(d.x, d.y), vt = (d.x * s->m->vel.y - d.y * s->m->vel.x) / r;
  double cth = -d.y / r;
  (void)w;

  return 0.5 * s->mass * vt * vt + s->mass * s->g * r * (1.0 - cth);
}

static double e_total(const pil_world *w, void *user) {
  return pil_world_energy(w) - ((swing_ref *)user)->E0;
}

static double e_rest(const pil_world *w, void *user) {
  return e_total(w, user) - e_swing(w, user);
}

int main(int argc, char **argv) {
  // data of problem
  const double g = 9.81;
  const double m = 1.0;
  const double M = 2.0;
  const double rho = 0.032;
  const double I = 2.56e-5;
  const double rhoR = 0.034;
  const double IR = 2.31e-5;
  const double n = 1.22;
  const double mc = 0.05;
  const double rstar = 0.4;
  const double th0 = -0.55;

  // Effective inertia of everything rope drags and spring that makes bounce
  // twice as fast as swing: k = 4 A g / r*
  const double A = m + M + 2 * I / (rho * rho) + 2 * IR * n * n / (rho * rho) +
                   mc * (n * rhoR / rho) * (n * rhoR / rho);
  const double k = 4 * A * g / rstar;

  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, -g);

  // Layout (metres)
  const double yt = 0.9;
  const double x1 = 1.05;
  const double x2 = 2.45;
  const double ycap = 0.46;
  const double xc0 = 1.75;
  pil_add_ground(w, (pil_vec2){x2 + 0.5, yt + 0.05},
                 (pil_vec2){0.4, yt + 0.05});
  pil_add_ground(w, (pil_vec2){0.4, -0.4}, (pil_vec2){x2 + 0.5, -0.4});

  // Pulleys: discs on fixed pivots, their tops at rope height
  pil_vec2 O1 = {x1, yt - rho};
  pil_vec2 O2 = {x2, yt - rho};
  pil_body *P1 = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(rho),
                                                  .mass = 0.05,
                                                  .inertia = I,
                                                  .pos = O1,
                                                  .flags = PIL_HIDE_COM,
                                                  .label = "P_1",
                                                  .label_dx = 20,
                                                  .label_dy = -22});
  pil_body *P2 = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(rho),
                                                  .mass = 0.05,
                                                  .inertia = I,
                                                  .pos = O2,
                                                  .flags = PIL_HIDE_COM,
                                                  .label = "P_2",
                                                  .label_dx = 16,
                                                  .label_dy = 22});
  pil_add_revolute(w, pil_on(P1, 0, 0), pil_at(O1.x, O1.y));
  pil_add_revolute(w, pil_on(P2, 0, 0), pil_at(O2.x, O2.y));

  // Block M on spring
  const double yM = 0.2;
  const double hM = 0.13;
  pil_body *Mb =
      pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.34, 2 * hM),
                                       .mass = M,
                                       .pos = {x2 + rho, yM},
                                       .label = "M",
                                       .flags = PIL_HIDE_COM});
  pil_add_prismatic(w, pil_on(Mb, 0, 0), pil_at(x2 + rho, 0), (pil_vec2){0, 1});
  double L0 = (yM - hM) - (-0.4); // spring length at release
  double rest = L0 + (M - m) * g / k;
  pil_element *sp =
      pil_add_spring(w, pil_at(x2 + rho, -0.4), pil_on(Mb, 0, -hM), k, rest);
  pil_element_set_label_at(sp, "k", -24, 0);

  // Swinging ball, released from rest at \theta_0 with r = r*
  pil_vec2 pivot = {x1, yt};
  pil_vec2 mp = {pivot.x + rstar * sin(th0), pivot.y - rstar * cos(th0)};
  pil_body *ball =
      pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.05),
                                       .mass = m,
                                       .pos = mp,
                                       .label = "m",
                                       .label_dx = -16,
                                       .label_dy = -10});
  pil_vec2 top = {x2 + rho, yt - rho};
  double rope = rstar + (top.y - (yM + hM));
  pil_add_pulley_rope(w, pil_on(ball, 0, 0), pil_on(Mb, 0, hM), pivot, top,
                      rope);

  // Inertia carried by rope: pulleys turn with rope speed
  pil_add_link(w, P1, PIL_COORD_ANGLE, Mb, PIL_COORD_Y, 1.0 / rho);
  pil_add_link(w, P2, PIL_COORD_ANGLE, Mb, PIL_COORD_Y, 1.0 / rho);

  // Belt -> pinch rollers -> capsule
  pil_vec2 R1c = {xc0, ycap + 0.025 + rhoR};
  pil_vec2 R2c = {xc0, ycap - 0.025 - rhoR};
  pil_body *R1 =
      pil_add_body(w, &(pil_body_desc){.shape = pil_disc(rhoR),
                                       .mass = 0.05,
                                       .inertia = IR,
                                       .pos = R1c,
                                       .flags = PIL_HIDE_COM | PIL_SPOKE});
  pil_body *R2 =
      pil_add_body(w, &(pil_body_desc){.shape = pil_disc(rhoR),
                                       .mass = 0.05,
                                       .inertia = IR,
                                       .pos = R2c,
                                       .flags = PIL_HIDE_COM | PIL_SPOKE});
  pil_add_revolute(w, pil_on(R1, 0, 0), pil_at(R1c.x, R1c.y));
  pil_add_revolute(w, pil_on(R2, 0, 0), pil_at(R2c.x, R2c.y));
  pil_add_link(w, R1, PIL_COORD_ANGLE, P2, PIL_COORD_ANGLE, n);
  pil_add_link(w, R2, PIL_COORD_ANGLE, P2, PIL_COORD_ANGLE, -n);

  pil_body *cap =
      pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.16, 0.05),
                                       .mass = mc,
                                       .pos = {xc0, ycap},
                                       .label = "m_c",
                                       .label_dx = 40,
                                       .label_dy = -6,
                                       .flags = PIL_HIDE_COM});
  pil_add_prismatic(w, pil_on(cap, 0, 0), pil_at(0, ycap), (pil_vec2){1, 0});
  pil_add_link(w, cap, PIL_COORD_X, R1, PIL_COORD_ANGLE, rhoR);
  pil_add_line(w, pil_on(P2, 0, 0), pil_on(R1, 0, 0), 0); // belt

  // Annotations
  pil_add_line(w, pil_at(pivot.x, pivot.y), pil_at(pivot.x, yM - 0.18), 1);
  pil_add_angle_between(w, pil_at(pivot.x, pivot.y), pil_at(pivot.x, 0),
                        pil_on(ball, 0, 0), 56, "\xcf\x91");
  pil_add_dimension(w, pil_at(pivot.x, pivot.y), pil_on(ball, 0, 0), -26, "r");
  pil_add_line(w, pil_at(xc0 - 0.3, ycap - 0.34), pil_at(xc0 - 0.3, ycap), 1);
  pil_add_dimension(w, pil_at(xc0 - 0.3, ycap - 0.28), pil_on(cap, 0, -0.28), 0,
                    "x_c");
  pil_add_gravity_indicator(w, (pil_vec2){0.5, 0.9});

  // Energy panel
  swing_ref ref = {ball, pivot, m, g, 0.0};
  ref.E0 = pil_world_energy(w) - e_swing(w, &ref);
  pil_plot *p = pil_add_plot(w, 0.0, 0.62, 1.0, 1.0, "t [s]", "E [J]", 52.0);
  pil_plot_add_signal(p, "E", e_total, &ref);
  pil_plot_add_signal(p, "E_s  swing", e_swing, &ref);
  pil_plot_add_signal(p, "E - E_s  bounce + shuttle", e_rest, &ref);
  pil_plots_record(w, 0.1, 5e-4);

  pil_view view;
  pil_view_fit_region(&view, 0.3, -0.5, x2 + 0.7, yt + 0.2, 900, 780, 0.0, 0.0,
                      1.0, 0.62, 18);
  int rc = example_main(w, &view, "Swinging Atwood machine with a shuttle",
                        argc, argv, "swing.png");

  pil_world_free(w);

  return rc;
}
