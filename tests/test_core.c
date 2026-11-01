/* Numerical regression tests against analytic solutions */

#include <math.h>
#include <pil/pil_core.h>
#include <stdio.h>

static int fails = 0;
static void check(const char *name, double got, double want, double tol) {
  int ok = fabs(got - want) <= tol;
  printf("  [%s] %-46s got % .6f  want % .6f  (tol %.1e)\n", ok ? "ok" : "FAIL",
         name, got, want, tol);
  if (!ok) {
    fails++;
  }
}

// 1. mass on a spring: period 2*pi*sqrt(m/k), energy conserved
static void test_spring(void) {
  puts("spring-mass");
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  pil_body *m = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.2, 0.2),
                                                 .mass = 2.0,
                                                 .pos = {1.0, 0}});
  pil_add_prismatic(w, pil_on(m, 0, 0), pil_at(0, 0), (pil_vec2){1, 0});
  pil_add_spring(w, pil_at(0, 0), pil_on(m, 0, 0), 50.0, 0.0);
  double T = 2 * M_PI * sqrt(2.0 / 50.0), E0 = pil_world_energy(w);
  pil_world_advance(w, T, 1e-3);
  check("x after one period", m->pos.x, 1.0, 1e-6);
  check("energy drift", pil_world_energy(w) / E0, 1.0, 1e-8);
  pil_world_free(w);
}

// 2. rigid rod pendulum: small-angle period of a physical pendulum, energy
static void test_pendulum(void) {
  puts("physical pendulum (revolute)");
  pil_world *w = pil_world_new();
  double L = 1.0, m = 3.0, th0 = 0.05;
  // uniform bar hinged at its end; COM at L/2
  pil_body *bar = pil_add_body(
      w, &(pil_body_desc){.shape = pil_rect(L, 0.05),
                          .mass = m,
                          .pos = {0.5 * L * sin(th0), -0.5 * L * cos(th0)},
                          .angle = th0 - M_PI / 2});
  pil_add_revolute(w, pil_at(0, 0), pil_on(bar, -0.5 * L, 0));
  double I = m * L * L / 3.0 + 0.0;    // bar about the hinge (thin)
  I = bar->inertia + m * 0.25 * L * L; // exact, includes thickness
  double T = 2 * M_PI * sqrt(I / (m * 9.81 * 0.5 * L));
  double E0 = pil_world_energy(w);
  pil_world_advance(w, T, 1e-3);
  check("angle after one period", bar->angle + M_PI / 2, th0, 2e-4);
  check("energy drift", pil_world_energy(w) / E0, 1.0, 1e-3);
  check("constraint error", pil_world_constraint_error(w), 0.0, 1e-6);
  pil_world_free(w);
}

// 3. cylinder rolling down an incline: a = g \sin(a) / 1.5
static void test_rolling(void) {
  puts("cylinder rolling on incline (rolling constraint)");
  pil_world *w = pil_world_new();
  double al = 0.4, R = 0.2;
  pil_vec2 u = {cos(al), sin(al)};
  pil_vec2 n = {-sin(al), cos(al)};
  pil_vec2 c0 = {R * n.x, R * n.y};
  pil_body *c = pil_add_body(
      w, &(pil_body_desc){.shape = pil_disc(R), .mass = 1.0, .pos = c0});
  pil_add_rolling(w, c, (pil_vec2){0, 0}, R, (pil_vec2){0, 0}, u);
  double t = 1.0;
  pil_world_advance(w, t, 1e-3);
  double a = 9.81 * sin(al) / 1.5;
  double s = c->pos.x * u.x + c->pos.y * u.y;
  check("distance along slope", s - (c0.x * u.x + c0.y * u.y), -0.5 * a * t * t,
        1e-6);
  check("no slip: v + R*\\omega",
        (c->vel.x * u.x + c->vel.y * u.y) + R * c->omega, 0.0, 1e-6);
  check("stays on plane", c->pos.x * n.x + c->pos.y * n.y, R, 1e-7);
  pil_world_free(w);
}

// 4. rope: taut -> acts like a rod; goes slack when pushed together
static void test_rope(void) {
  puts("rope");
  pil_world *w = pil_world_new();
  pil_body *b = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.05),
                                                 .mass = 1.0,
                                                 .pos = {1, 0}});
  pil_element *r = pil_add_distance(w, pil_at(0, 0), pil_on(b, 0, 0),
                                    pil_LENGTH_INITIAL, pil_ROPE);
  pil_world_advance(w, 0.001, 1e-3);
  // hanging sideways, released from rest: pendulum, tension at start = 0 +
  // centripetal 0 -> m g \sin...
  check("initial tension (horizontal rope, v=0)", pil_element_force(r), 0.0,
        1e-3);
  pil_world_advance(w, 0.5, 1e-3);
  double len = sqrt(b->pos.x * b->pos.x + b->pos.y * b->pos.y);
  check("length preserved while taut", len, 1.0, 1e-5);
  pil_world_free(w);

  w = pil_world_new();
  b = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.05),
                                       .mass = 1.0,
                                       .pos = {0.5, 0}});
  pil_add_distance(w, pil_at(0, 0), pil_on(b, 0, 0), 1.0,
                   pil_ROPE); // slack: 0.5 < 1.0
  pil_world_advance(w, 0.2, 1e-3);
  check("slack rope: free fall y = -g t^2/2", b->pos.y, -0.5 * 9.81 * 0.04,
        1e-5);
  pil_world_free(w);
}

// 5. Atwood-like: two masses on a rope over a fixed point via two ropes
//    use inclined slot + spring + friction-free block: a = g \sin(al)
static void test_slot(void) {
  puts("block on frictionless slot (prismatic)");
  pil_world *w = pil_world_new();
  double al = 0.6;
  pil_vec2 u = {cos(al), sin(al)};
  pil_body *b = pil_add_body(
      w,
      &(pil_body_desc){.shape = pil_rect(0.2, 0.1), .mass = 2.0, .angle = al});
  pil_add_prismatic(w, pil_on(b, 0, 0), pil_at(0, 0), u);
  pil_world_advance(w, 0.8, 1e-3);
  double s = b->pos.x * u.x + b->pos.y * u.y;
  check("s = -g sin(a) t^2 / 2", s, -0.5 * 9.81 * sin(al) * 0.64, 1e-6);
  check("angle held", b->angle, al, 1e-7);
  pil_world_free(w);
}

// 6. block resting on a floor (penalty contact): settles with N = m g
static void test_floor(void) {
  puts("block on floor (contact + friction)");
  pil_world *w = pil_world_new();
  pil_body *b = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.4, 0.2),
                                                 .mass = 5.0,
                                                 .pos = {0, 0.1},
                                                 .vel = {1.0, 0}});
  pil_add_floor_contact(w, b, pil_at(0, 0), (pil_vec2){1, 0}, 2e5, 2e3, 0.3);
  pil_world_advance(w, 0.2, 1e-4);
  // sliding friction decelerates at \mu g = 2.943 m/s^2
  check("sliding deceleration", 1.0 - b->vel.x, 0.3 * 9.81 * 0.2, 2e-3);
  check("rest height (small penetration)", b->pos.y, 0.1, 2e-3);
  pil_world_free(w);
}

// 7. double pendulum chain: constraint drift + energy over 10 s
static void test_chain(void) {
  puts("double pendulum (chaotic regime, 10 s)");
  pil_world *w = pil_world_new();
  pil_body *a = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(1.0, 0.04),
                                                 .mass = 1.0,
                                                 .pos = {0.5, 0}});
  pil_body *b = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(1.0, 0.04),
                                                 .mass = 1.0,
                                                 .pos = {1.5, 0}});
  pil_add_revolute(w, pil_at(0, 0), pil_on(a, -0.5, 0));
  pil_add_revolute(w, pil_on(a, 0.5, 0), pil_on(b, -0.5, 0));
  double E0 = pil_world_energy(w);
  pil_world_advance(w, 10.0, 5e-4);
  check("energy drift / (m g L)", (pil_world_energy(w) - E0) / (2 * 9.81), 0.0,
        5e-4);
  check("constraint error", pil_world_constraint_error(w), 0.0, 1e-5);
  pil_world_free(w);
}

// 8. Atwood machine: a = g (m1-m2)/(m1+m2), T = 2 m1 m2 g /(m1+m2); rope goes
// slack if both are lifted
static void test_atwood(void) {
  puts("Atwood machine (pulley rope)");
  pil_world *w = pil_world_new();
  double m1 = 3.0, m2 = 1.0, R = 0.15;
  pil_body *A = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.2, 0.2),
                                                 .mass = m1,
                                                 .pos = {-R, 1.0}});
  pil_body *B = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.2, 0.2),
                                                 .mass = m2,
                                                 .pos = {+R, 1.0}});
  pil_add_prismatic(w, pil_on(A, 0, 0), pil_at(-R, 0), (pil_vec2){0, 1});
  pil_add_prismatic(w, pil_on(B, 0, 0), pil_at(+R, 0), (pil_vec2){0, 1});
  pil_element *rope = pil_add_pulley_rope(
      w, pil_on(A, 0, 0.1), pil_on(B, 0, 0.1), (pil_vec2){-R, 2.5},
      (pil_vec2){+R, 2.5}, pil_LENGTH_INITIAL);
  double a = 9.81 * (m1 - m2) / (m1 + m2), t = 0.5;
  pil_world_advance(w, t, 1e-3);
  check("m1 falls 1/2 a t^2", A->pos.y - 1.0, -0.5 * a * t * t, 1e-6);
  check("m2 rises 1/2 a t^2", B->pos.y - 1.0, +0.5 * a * t * t, 1e-6);
  check("tension", pil_element_force(rope), 2 * m1 * m2 * 9.81 / (m1 + m2),
        1e-3);
  check("constraint error", pil_world_constraint_error(w), 0.0, 1e-8);
  pil_world_free(w);

  // equal masses at rest: tension = m g, nothing moves
  w = pil_world_new();
  A = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.2, 0.2),
                                       .mass = 2.0,
                                       .pos = {-R, 1.0}});
  B = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.2, 0.2),
                                       .mass = 2.0,
                                       .pos = {+R, 1.0}});
  pil_add_prismatic(w, pil_on(A, 0, 0), pil_at(-R, 0), (pil_vec2){0, 1});
  pil_add_prismatic(w, pil_on(B, 0, 0), pil_at(+R, 0), (pil_vec2){0, 1});
  rope = pil_add_pulley_rope(w, pil_on(A, 0, 0.1), pil_on(B, 0, 0.1),
                             (pil_vec2){-R, 2.5}, (pil_vec2){+R, 2.5},
                             pil_LENGTH_INITIAL);
  pil_world_advance(w, 1.0, 1e-3);
  check("balanced: no motion", A->pos.y, 1.0, 1e-7);
  check("balanced: tension = m g", pil_element_force(rope), 2.0 * 9.81, 1e-3);
  pil_world_free(w);
}

// 9. mouse grab: pick a body, drag its target; the grabbed point follows
static void test_grab(void) {
  puts("mouse grab (pick + soft spring)");
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  pil_body *b = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.4, 0.2),
                                                 .mass = 2.0,
                                                 .pos = {0, 0}});
  pil_vec2 local;
  check("pick finds the body",
        pil_world_pick(w, (pil_vec2){0.1, 0.05}, &local) == b, 1.0, 0);
  check("pick misses empty space",
        pil_world_pick(w, (pil_vec2){0.5, 0.5}, NULL) == NULL, 1.0, 0);

  // grab at the COM: pure translation, critically damped -> COM lands on the
  // target
  pil_world_set_grab(w, b, (pil_vec2){0, 0}, (pil_vec2){1.0, 0.5},
                     300.0 * b->mass, 2.0 * sqrt(300.0) * b->mass);
  pil_world_advance(w, 3.0, 1e-3);
  check("COM grab: x at target", b->pos.x, 1.0, 1e-3);
  check("COM grab: y at target", b->pos.y, 0.5, 1e-3);
  check("COM grab: no rotation", b->angle, 0.0, 1e-9);
  pil_world_clear_grab(w);

  // off-centre grab: the body rotates, but the grabbed point reaches the target
  pil_vec2 loc = {0.15, 0.05}, target = {2.0, 1.0};
  pil_world_set_grab(w, b, loc, target, 300.0 * b->mass,
                     2.0 * sqrt(300.0) * b->mass);
  pil_world_advance(w, 8.0, 1e-3);
  pil_vec2 p = pil_body_point(b, loc);
  check("off-centre grab: point x at target", p.x, target.x, 5e-3);
  check("off-centre grab: point y at target", p.y, target.y, 5e-3);
  pil_world_clear_grab(w);
  pil_world_free(w);
}

int main(void) {
  test_spring();
  test_pendulum();
  test_rolling();
  test_rope();
  test_slot();
  test_floor();
  test_chain();
  test_atwood();
  test_grab();
  printf("\n%s (%d failure%s)\n", fails ? "FAILED" : "ALL PASSED", fails,
         fails == 1 ? "" : "s");
  return fails ? 1 : 0;
}
