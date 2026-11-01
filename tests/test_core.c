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
  double T = 2 * M_PI * sqrt(2.0 / 50.0);
  double E0 = pil_world_energy(w);
  pil_world_advance(w, T, 1e-3);
  check("x after one period", m->pos.x, 1.0, 1e-6);
  check("energy drift", pil_world_energy(w) / E0, 1.0, 1e-8);
  pil_world_free(w);
}

// 2. rigid rod pendulum: small-angle period of a physical pendulum, energy
static void test_pendulum(void) {
  puts("physical pendulum (revolute)");
  pil_world *w = pil_world_new();
  double L = 1.0;
  double m = 3.0;
  double th0 = 0.05;
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
  double al = 0.4;
  double R = 0.2;
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
                                    PIL_LENGTH_INITIAL, PIL_ROPE);
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
                   PIL_ROPE); // slack: 0.5 < 1.0
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

// 8. Atwood machine:
// a = g (m1-m2)/(m1+m2)
// T = 2 m1 m2 g /(m1+m2)
// rope goes slack if both are lifted
static void test_atwood(void) {
  puts("Atwood machine (pulley rope)");
  pil_world *w = pil_world_new();
  double m1 = 3.0;
  double m2 = 1.0;
  double R = 0.15;
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
      (pil_vec2){+R, 2.5}, PIL_LENGTH_INITIAL);
  double a = 9.81 * (m1 - m2) / (m1 + m2);
  double t = 0.5;
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
                             PIL_LENGTH_INITIAL);
  pil_world_advance(w, 1.0, 1e-3);
  check("balanced: no motion", A->pos.y, 1.0, 1e-7);
  check("balanced: tension = m g", pil_element_force(rope), 2.0 * 9.81, 1e-3);
  pil_world_free(w);
}

// 9. mouse grab: pick a body, drag its target; grabbed point follows
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

  // grab at COM: pure translation, critically damped -> COM lands on target
  pil_world_set_grab(w, b, (pil_vec2){0, 0}, (pil_vec2){1.0, 0.5},
                     300.0 * b->mass, 2.0 * sqrt(300.0) * b->mass);
  pil_world_advance(w, 3.0, 1e-3);
  check("COM grab: x at target", b->pos.x, 1.0, 1e-3);
  check("COM grab: y at target", b->pos.y, 0.5, 1e-3);
  check("COM grab: no rotation", b->angle, 0.0, 1e-9);
  pil_world_clear_grab(w);

  // off-centre grab: the body rotates, but the grabbed point reaches the target
  pil_vec2 loc = {0.15, 0.05};
  pil_vec2 target = {2.0, 1.0};
  pil_world_set_grab(w, b, loc, target, 300.0 * b->mass,
                     2.0 * sqrt(300.0) * b->mass);
  pil_world_advance(w, 8.0, 1e-3);
  pil_vec2 p = pil_body_point(b, loc);
  check("off-centre grab: point x at target", p.x, target.x, 5e-3);
  check("off-centre grab: point y at target", p.y, target.y, 5e-3);
  pil_world_clear_grab(w);
  pil_world_free(w);
}

// 10. circle rolling inside a fixed bowl:
// T = 2 \pi \sqrt((1 + I/(m r^2)) (R - r) / g)
static void test_ball_in_bowl(void) {
  puts("ball rolling inside a bowl (circle rolling, fixed track)");
  pil_world *w = pil_world_new();
  double R = 0.5;
  double r = 0.1;
  double D = R - r;
  double psi0 = -M_PI / 2 + 0.02;
  pil_vec2 c = {0, 0.6};
  pil_body *ball = pil_add_body(
      w, &(pil_body_desc){.shape = pil_disc(r),
                          .mass = 1.0,
                          .pos = {c.x + D * cos(psi0), c.y + D * sin(psi0)}});
  pil_add_circle_rolling(w, pil_on(ball, 0, 0), r, pil_at(c.x, c.y), R, 1);
  double T = 2 * M_PI * sqrt(1.5 * D / 9.81);
  double E0 = pil_world_energy(w);
  pil_vec2 p0 = ball->pos;
  pil_world_advance(w, T, 1e-3);
  check("returns after one analytic period (x)", ball->pos.x, p0.x, 2e-5);
  check("energy drift / |E|", (pil_world_energy(w) - E0) / fabs(E0), 0.0, 1e-7);
  check("constraint error", pil_world_constraint_error(w), 0.0, 1e-6);
  // no slip: contact-point velocity along track is zero => \omega = -(D
  // psidot)/r
  pil_world_advance(w, 0.3, 1e-3);
  pil_vec2 d = {ball->pos.x - c.x, ball->pos.y - c.y};
  double psid =
      (d.x * ball->vel.y - d.y * ball->vel.x) / (d.x * d.x + d.y * d.y);
  check("no slip: r*omega + D*psidot", r * ball->omega + D * psid, 0.0, 1e-6);
  pil_world_free(w);
}

// 11. ball in a MOVING bowl: horizontal momentum and energy are conserved
static void test_ball_in_moving_bowl(void) {
  puts("ball in a movable bowl (both bodies dynamic)");
  pil_world *w = pil_world_new();
  double R = 0.5;
  double r = 0.1;
  double D = R - r;
  double psi0 = -M_PI / 2 + 0.6;
  double mb = 2.0;
  double mr = 1.0;
  pil_vec2 c = {0, 0.6};
  pil_body *bowl = pil_add_body(
      w, &(pil_body_desc){
             .shape = pil_disc(R), .mass = mb, .pos = c, .flags = PIL_NO_FILL});
  pil_body *ball = pil_add_body(
      w, &(pil_body_desc){.shape = pil_disc(r),
                          .mass = mr,
                          .pos = {c.x + D * cos(psi0), c.y + D * sin(psi0)}});
  pil_add_prismatic(w, pil_on(bowl, 0, 0), pil_at(0, c.y),
                    (pil_vec2){1, 0}); /* bowl slides horizontally */
  pil_add_fixed_angle(w, bowl, NULL, 0.0);
  pil_add_circle_rolling(w, pil_on(ball, 0, 0), r, pil_on(bowl, 0, 0), R, 1);
  double E0 = pil_world_energy(w);
  pil_world_advance(w, 2.0, 5e-4);
  check("horizontal momentum", mb * bowl->vel.x + mr * ball->vel.x, 0.0, 1e-6);
  check("energy drift / |E|", (pil_world_energy(w) - E0) / fabs(E0), 0.0, 1e-5);
  check("bowl actually moved", fabs(bowl->pos.x) > 0.01, 1.0, 0);
  check("constraint error", pil_world_constraint_error(w), 0.0, 1e-5);
  pil_world_free(w);
}

// 12. ball over top of a cylinder (convex contact)
static void test_ball_over_cylinder(void) {
  puts("ball on top of a cylinder (convex contact)");
  pil_world *w = pil_world_new();
  double R = 0.5;
  double r = 0.1;
  double D = R + r;
  double psi0 = M_PI / 2 - 0.1;
  pil_body *ball =
      pil_add_body(w, &(pil_body_desc){.shape = pil_disc(r),
                                       .mass = 1.0,
                                       .pos = {D * cos(psi0), D * sin(psi0)}});
  pil_add_circle_rolling(w, pil_on(ball, 0, 0), r, pil_at(0, 0), R, 0);
  double E0 = pil_world_energy(w);
  pil_world_advance(w, 0.5, 5e-4);
  check("energy drift / |E|", (pil_world_energy(w) - E0) / fabs(E0), 0.0, 1e-6);
  check("constraint error", pil_world_constraint_error(w), 0.0, 1e-6);
  check("rolled down", ball->pos.x > 0.0 + D * cos(psi0) + 0.05, 1.0, 0);
  pil_world_free(w);
}

// 13. weld: compound pendulum (bar + welded end mass), T from I_total/(m g d)
static void test_weld(void) {
  puts("weld (compound pendulum)");
  pil_world *w = pil_world_new();
  double L = 1.0;
  double m = 2.0;
  double M = 1.5;
  double th0 = 0.04;
  pil_vec2 dir = {sin(th0), -cos(th0)};
  pil_body *bar = pil_add_body(
      w, &(pil_body_desc){.shape = pil_rect(L, 0.04),
                          .mass = m,
                          .pos = {0.5 * L * dir.x, 0.5 * L * dir.y},
                          .angle = th0 - M_PI / 2});
  pil_body *bob =
      pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.05),
                                       .mass = M,
                                       .pos = {L * dir.x, L * dir.y},
                                       .angle = th0 - M_PI / 2});
  pil_add_revolute(w, pil_at(0, 0), pil_on(bar, -0.5 * L, 0));
  pil_add_weld(w, pil_on(bar, 0.5 * L, 0), pil_on(bob, 0, 0));
  double I = bar->inertia + m * 0.25 * L * L + bob->inertia + M * L * L;
  double T = 2 * M_PI * sqrt(I / (9.81 * (m * 0.5 * L + M * L)));
  pil_world_advance(w, T, 5e-4);
  check("compound period (angle returns)", bar->angle + M_PI / 2, th0, 3e-4);
  check("bob stays welded (distance)",
        hypot(bob->pos.x - pil_body_point(bar, (pil_vec2){0.5 * L, 0}).x,
              bob->pos.y - pil_body_point(bar, (pil_vec2){0.5 * L, 0}).y),
        0.0, 1e-6);
  check("same angle as the bar", bob->angle - bar->angle, 0.0, 1e-6);
  pil_world_free(w);
}

// 14. rope snapping taut at speed: no visible stretch and snap is nelastic
// (energy never rises)
static void test_rope_snap(void) {
  puts("rope snapping taut at ~5 m/s");
  pil_world *w = pil_world_new();
  pil_body *b = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.05),
                                                 .mass = 1.0,
                                                 .pos = {0.3, 0.5}});
  pil_add_distance(w, pil_at(0, 0), pil_on(b, 0, 0), 1.0,
                   PIL_ROPE); // slack: starts 0.58 m away
  double E0 = pil_world_energy(w), worst = 0, Emax = E0;
  for (int i = 0; i < 1500; i++) {
    pil_world_step(w, 5e-4);
    double e = pil_world_constraint_error(w), E = pil_world_energy(w);
    if (e > worst) {
      worst = e;
    }
    if (E > Emax) {
      Emax = E;
    }
  }
  check("max rope stretch", worst, 0.0, 2e-3);
  check("energy never rises (inelastic snap)", Emax - E0, 0.0, 1e-6);
  check("energy was dissipated by the snap", pil_world_energy(w) < E0 - 0.1,
        1.0, 0);
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
  test_ball_in_bowl();
  test_ball_in_moving_bowl();
  test_ball_over_cylinder();
  test_weld();
  test_rope_snap();
  printf("\n%s (%d failure%s)\n", fails ? "FAILED" : "ALL PASSED", fails,
         fails == 1 ? "" : "s");
  return fails ? 1 : 0;
}
