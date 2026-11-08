/* Numerical regression tests against analytic solutions */

#include <math.h>
#include <pil/pil_circuit.h>
#include <pil/pil_core.h>
#include <pil/pil_lagrange.h>
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
  double I = m * L * L / 3.0 + 0.0;    // bar about hinge (thin)
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

static double sig_y(const pil_world *w, void *user) {
  (void)w;

  return ((pil_body *)user)->pos.y;
}

static void test_link_and_plot(void) {
  puts("linear link (effective inertia) and plot recording");
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, -9.81);
  double ma = 1.0;
  double mb = 0.5;
  double ratio = -2.0;
  pil_body *a = pil_add_body(
      w,
      &(pil_body_desc){.shape = pil_rect(0.2, 0.2), .mass = ma, .pos = {0, 5}});
  pil_body *b = pil_add_body(
      w,
      &(pil_body_desc){.shape = pil_rect(0.2, 0.2), .mass = mb, .pos = {1, 0}});
  pil_add_prismatic(w, pil_on(a, 0, 0), pil_at(0, 0), (pil_vec2){0, 1});
  pil_add_prismatic(w, pil_on(b, 0, 0), pil_at(0, 0), (pil_vec2){1, 0});
  // b.x = ratio * a.y + const: b moves |ratio| times as far as a drops
  pil_add_link(w, b, PIL_COORD_X, a, PIL_COORD_Y, ratio);
  pil_plot *p = pil_add_plot(w, 0.5, 0.5, 1, 1, "t [s]", "y", 1.0);
  pil_plot_add_signal(p, "y", sig_y, a);
  pil_plots_record(w, 0.1, 1e-3);
  check("recording leaves the state untouched", a->pos.y, 5.0, 1e-12);
  check("recording rewinds the clock", pil_world_time(w), 0.0, 1e-12);
  // Effective mass ma + ratio^2 mb: a_y = -g ma / (ma + ratio^2 mb)
  double acc = 9.81 * ma / (ma + ratio * ratio * mb);
  check("plot sample at t = 1 s (y = 5 - acc/2)", pil_plot_value(p, 0, 10),
        5.0 - 0.5 * acc, 2e-4);
  pil_world_advance(w, 1.0, 5e-4);
  check("fall with effective inertia", a->pos.y, 5.0 - 0.5 * acc, 2e-4);
  check("follower moves 2x", b->pos.x - 1.0, ratio * (a->pos.y - 5.0), 1e-6);
  pil_world_free(w);
}

static void test_circuits(void) {
  puts("circuits (RC, RL, RLC, LC energy, AC, diode, switch)");
  // RC charge: v = V (1 - \exp(-t/RC))
  {
    pil_world *w = pil_world_new();
    pil_circuit *c = pil_add_circuit(w);
    pil_node a = pil_circuit_node(c, 0, 1), b = pil_circuit_node(c, 1, 1),
             g = pil_circuit_node(c, 0, 0);
    pil_circuit_ground(c, g);
    pil_vsource_dc(c, a, g, 5.0, "V");
    pil_resistor(c, a, b, 1e3, "R");
    pil_celement *cap = pil_capacitor(c, b, g, 1e-6, 0.0, "C");
    pil_world_advance(w, 2e-3, 1e-3);
    check("RC charge at 2 tau", pil_celement_voltage(cap),
          5.0 * (1 - exp(-2.0)), 5e-4);
    pil_world_free(w);
  }
  // RL: i = V/R (1 - \exp(-tR/L))
  {
    pil_world *w = pil_world_new();
    pil_circuit *c = pil_add_circuit(w);
    pil_node a = pil_circuit_node(c, 0, 1), b = pil_circuit_node(c, 1, 1),
             g = pil_circuit_node(c, 0, 0);
    pil_circuit_ground(c, g);
    pil_vsource_dc(c, a, g, 10.0, "V");
    pil_resistor(c, a, b, 10.0, "R");
    pil_celement *ind = pil_inductor(c, b, g, 0.1, 0.0, "L");
    pil_world_advance(w, 20e-3, 1e-3);
    check("RL current at 2 tau", pil_celement_current(ind),
          1.0 * (1 - exp(-2.0)), 5e-4);
    pil_world_free(w);
  }
  // series RLC, capacitor initially charged
  {
    pil_world *w = pil_world_new();
    pil_circuit *c = pil_add_circuit(w);
    pil_circuit_set_step(c, 1e-6);
    pil_node a = pil_circuit_node(c, 0, 1), b = pil_circuit_node(c, 1, 1),
             g = pil_circuit_node(c, 0, 0);
    pil_circuit_ground(c, g);
    double R = 2.0;
    double L = 1e-3;
    double C = 1e-6;
    pil_celement *cap = pil_capacitor(c, a, g, C, 1.0, "C");
    pil_inductor(c, a, b, L, 0.0, "L");
    pil_resistor(c, b, g, R, "R");
    double al = R / (2 * L);
    double w0 = 1 / sqrt(L * C);
    double wd = sqrt(w0 * w0 - al * al);
    double t = 150e-6;
    pil_world_advance(w, t, 1e-5);
    check("RLC ringing, v_C(150 us)", pil_celement_voltage(cap),
          exp(-al * t) * (cos(wd * t) + al / wd * sin(wd * t)), 2e-3);
    pil_world_free(w);
  }
  // LC tank: trapezoidal rule keeps energy
  {
    pil_world *w = pil_world_new();
    pil_circuit *c = pil_add_circuit(w);
    pil_node a = pil_circuit_node(c, 0, 1), g = pil_circuit_node(c, 0, 0);
    pil_circuit_ground(c, g);
    pil_capacitor(c, a, g, 1e-6, 2.0, "C");
    pil_inductor(c, a, g, 1e-3, 0.0, "L");
    double E0 = pil_circuit_energy(c);
    pil_world_advance(w, 5e-3, 1e-4);
    check("LC energy after 30 cycles", pil_circuit_energy(c) / E0, 1.0, 1e-3);
    pil_world_free(w);
  }
  // RC low-pass: gain 1/\sqrt(2) at w = 1/RC
  {
    pil_world *w = pil_world_new();
    pil_circuit *c = pil_add_circuit(w);
    pil_node a = pil_circuit_node(c, 0, 1), b = pil_circuit_node(c, 1, 1),
             g = pil_circuit_node(c, 0, 0);
    pil_circuit_ground(c, g);
    pil_vsource_ac(c, a, g, 1.0, 1000.0 / (2 * M_PI), 0.0, "v");
    pil_resistor(c, a, b, 1e3, "R");
    pil_celement *cap = pil_capacitor(c, b, g, 1e-6, 0.0, "C");
    pil_world_advance(w, 50e-3, 1e-3);
    double peak = 0;
    for (int i = 0; i < 400; i++) {
      pil_world_advance(w, 2e-5, 1e-3);
      peak = fmax(peak, fabs(pil_celement_voltage(cap)));
    }

    check("RC low-pass gain at the corner", peak, 1 / sqrt(2.0), 2e-3);
    pil_world_free(w);
  }
  // half-wave rectifier: peak output = V_pk - V_d with V_d from diode law
  {
    pil_world *w = pil_world_new();
    pil_circuit *c = pil_add_circuit(w);
    pil_node a = pil_circuit_node(c, 0, 1), b = pil_circuit_node(c, 1, 1),
             g = pil_circuit_node(c, 0, 0);
    pil_circuit_ground(c, g);
    pil_vsource_ac(c, a, g, 5.0, 50.0, 0.0, "v");
    pil_diode(c, a, b, 1e-12, 0.026, "D");
    pil_resistor(c, b, g, 1e3, "R");
    double vd = 0.6;
    for (int i = 0; i < 60; i++) { // Newton on 5 = vd + R Is (\exp(vd/nvt) - 1)
      double f = vd + 1e3 * 1e-12 * (exp(vd / 0.026) - 1) - 5.0;
      double df = 1 + 1e3 * 1e-12 * exp(vd / 0.026) / 0.026;
      vd -= f / df;
    }

    double vmax = -1e9;
    double vmin = 1e9;
    int worst = 0;
    for (int i = 0; i < 400; i++) {
      pil_world_advance(w, 5e-5, 1e-3);
      double vo = pil_circuit_node_voltage(c, b);
      vmax = fmax(vmax, vo);
      vmin = fmin(vmin, vo);
      worst = pil_circuit_newton_iterations(c) > worst
                  ? pil_circuit_newton_iterations(c)
                  : worst;
    }

    check("rectifier peak = V_pk - V_d", vmax, 5.0 - vd, 3e-3);
    check("rectifier blocks the negative half", vmin, 0.0, 1e-6);
    check("diode Newton stays cheap (< 30 iterations)", worst < 30, 1, 0);
    pil_world_free(w);
  }
  // switch closes at 1 ms, then RC charges
  {
    pil_world *w = pil_world_new();
    pil_circuit *c = pil_add_circuit(w);
    pil_node a = pil_circuit_node(c, 0, 1), s = pil_circuit_node(c, 1, 1),
             b = pil_circuit_node(c, 2, 1), g = pil_circuit_node(c, 0, 0);
    pil_circuit_ground(c, g);
    pil_circuit_set_theta(c, 1.0); // backward Euler: no ringing at a switch
    pil_circuit_set_step(c, 2e-6);
    pil_vsource_dc(c, a, g, 5.0, "V");
    pil_switch(c, a, s, 1e-3, 0, "S");
    pil_resistor(c, s, b, 1e3, "R");
    pil_celement *cap = pil_capacitor(c, b, g, 1e-6, 0.0, "C");
    pil_world_advance(w, 0.9e-3, 1e-4);
    // open switch is 1 Gohm, so a few microvolts leak through
    check("open switch: no charge", pil_celement_voltage(cap), 0.0, 2e-5);
    pil_world_advance(w, 2e-3, 1e-4);
    check("one \\tau after closing", pil_celement_voltage(cap),
          5.0 * (1 - exp(-1.9)), 2e-2);
    pil_world_free(w);
  }
}

static void test_electromechanics(void) {
  puts("electromechanical coupling (DC motor, voice coil)");
  // motor spinning a flywheel against viscous drag:
  // \omega(t) = \omega_{ss} (1 - \exp(-t/\tau)),
  // \omega_{ss} = Kt V / (R c + Kt^2),
  // \tau = J R / (R c + Kt^2)
  {
    pil_world *w = pil_world_new();
    pil_world_set_gravity(w, 0, 0);
    double Kt = 0.05;
    double R = 1.0;
    double V = 5.0;
    double c = 1e-3;
    double J = 2e-4;
    pil_body *wheel = pil_add_body(
        w, &(pil_body_desc){.shape = pil_disc(0.1), .mass = 1, .inertia = J});
    pil_add_revolute(w, pil_on(wheel, 0, 0), pil_at(0, 0));
    pil_add_torsion_spring(w, wheel, NULL, 0.0, 0.0, c);
    pil_circuit *ck = pil_add_circuit(w);
    pil_node a = pil_circuit_node(ck, 0, 1), g = pil_circuit_node(ck, 0, 0);
    pil_circuit_ground(ck, g);
    pil_vsource_dc(ck, a, g, V, "V");
    pil_motor(ck, a, g, wheel, NULL, Kt, R, 1e-6, "M");
    double den = R * c + Kt * Kt, wss = Kt * V / den, tau = J * R / den;
    pil_world_advance(w, 0.2, 1e-4);
    check("flywheel speed at t = 0.2 s", wheel->omega,
          wss * (1 - exp(-0.2 / tau)), 0.01 * wss);
    pil_world_free(w);
  }
  // voice coil holding a spring-mounted mass: x_ss = K V / (R k)
  {
    pil_world *w = pil_world_new();
    pil_world_set_gravity(w, 0, 0);
    double K = 2.0;
    double R = 4.0;
    double V = 10.0;
    double k = 100.0;
    pil_body *m = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.1, 0.1),
                                                   .mass = 0.1,
                                                   .pos = {0, 0}});
    pil_add_prismatic(w, pil_on(m, 0, 0), pil_at(0, 0), (pil_vec2){1, 0});
    pil_add_spring(w, pil_at(-0.5, 0), pil_on(m, 0, 0), k, 0.5);
    pil_circuit *ck = pil_add_circuit(w);
    pil_node a = pil_circuit_node(ck, 0, 1), g = pil_circuit_node(ck, 0, 0);
    pil_circuit_ground(ck, g);
    pil_vsource_dc(ck, a, g, V, "V");
    pil_voice_coil(ck, a, g, m, NULL, (pil_vec2){1, 0}, K, R, 1e-5, "VC");
    pil_world_advance(w, 3.0, 2e-4);
    check("voice-coil displacement", m->pos.x, K * V / (R * k), 1e-4);
    pil_world_free(w);
  }
}

static void test_stops(void) {
  puts("impact stops");
  // free mass hits a wall at x = 1 with restitution e:
  // speed -> -e v and lost kinetic energy is (1 - e^2) m v^2 / 2
  double es[] = {1.0, 0.6, 0.0};
  for (int i = 0; i < 3; i++) {
    pil_world *w = pil_world_new();
    pil_world_set_gravity(w, 0, 0);
    double v0 = 2.0;
    double m = 1.5;
    pil_body *b = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.1, 0.1),
                                                   .mass = m,
                                                   .pos = {0, 0},
                                                   .vel = {v0, 0}});
    pil_add_prismatic(w, pil_on(b, 0, 0), pil_at(0, 0), (pil_vec2){1, 0});
    pil_add_stop(w, b, PIL_COORD_X, 1.0, 1, es[i]);
    pil_world_advance(w, 0.8, 1e-3);
    check("rebound speed", b->vel.x, -es[i] * v0, 1e-3);
    check("energy lost", pil_world_dissipated(w),
          0.5 * m * v0 * v0 * (1 - es[i] * es[i]), 1e-3);
    pil_world_free(w);
  }
}

static void test_circle_contact(void) {
  puts("circle contact (loop-the-loop)");
  double g = 9.81, Rl = 1.0, r = 0.05, D = Rl - r, m = 1.0;
  // frictionless disc sliding inside a loop leaves wall where \cos(\phi) = (2 g
  // D - v0^2) / (3 g D), \phi measured from bottom
  {
    pil_world *w = pil_world_new();
    pil_body *b =
        pil_add_body(w, &(pil_body_desc){.shape = pil_disc(r),
                                         .mass = m,
                                         .pos = {0, Rl - D},
                                         .vel = {sqrt(3 * g * D), 0}});
    pil_add_circle_contact(w, b, pil_at(0, Rl), Rl, 1, 4e4, 40.0, 0.0);
    double want = acos(-1.0 / 3.0), got = -1;
    for (int i = 0; i < 4000 && got < 0; i++) {
      pil_world_advance(w, 1e-3, 1e-4);
      pil_vec2 d = {b->pos.x, b->pos.y - Rl};
      if (D - hypot(d.x, d.y) > 2e-3) {
        got = atan2(d.x, -d.y);
      }
    }
    check("detachment angle", got, want, 0.12);
    pil_world_free(w);
  }
  // rolling without slipping (high mu): bottom speed^2 = 4/3 g h for a solid
  // disc released from rest at height h above bottom of track
  {
    pil_world *w = pil_world_new();
    double th = 1.0; // start angle from bottom
    pil_body *b = pil_add_body(
        w, &(pil_body_desc){.shape = pil_disc(r),
                            .mass = m,
                            .pos = {D * sin(th), Rl - D * cos(th)}});
    pil_add_circle_contact(w, b, pil_at(0, Rl), Rl, 1, 4e4, 40.0, 5.0);
    double h = D * (1 - cos(th)), vmax = 0;
    for (int i = 0; i < 3000; i++) {
      pil_world_advance(w, 1e-3, 1e-4);
      double sp = hypot(b->vel.x, b->vel.y);
      if (sp > vmax) {
        vmax = sp;
      }
    }
    check("rolling bottom speed", vmax, sqrt(4.0 / 3.0 * g * h), 0.03);
    pil_world_free(w);
  }
}

typedef struct {
  double m1, m2, l, k, L0, Om, g;
} tower_params;

// Block m1 (distance x along an arm turning at Om), pendulum m2 (length l,
// angle gamma, swinging across arm)
static double tower_T(const double *q, const double *qd, void *u) {
  const tower_params *p = u;
  double x = q[0], xd = qd[0], gd = qd[1], sg = sin(q[1]), cg = cos(q[1]);
  double t1 = 0.5 * p->m1 * (xd * xd + p->Om * p->Om * x * x);
  double vx = xd - p->Om * p->l * sg, vy = p->l * cg * gd + p->Om * x,
         vz = p->l * sg * gd;
  return t1 + 0.5 * p->m2 * (vx * vx + vy * vy + vz * vz);
}

static double tower_V(const double *q, void *u) {
  const tower_params *p = u;
  return 0.5 * p->k * (q[0] - p->L0) * (q[0] - p->L0) -
         p->m2 * p->g * p->l * cos(q[1]);
}

static double pend_T(const double *q, const double *qd, void *u) {
  (void)q;
  (void)u;
  return 0.5 * qd[0] * qd[0];
}

static double pend_V(const double *q, void *u) {
  (void)u;
  return -9.81 * cos(q[0]);
}

static void test_lagrangian(void) {
  puts("generalised coordinates (Euler-Lagrange)");
  // Pendulum, l = 1, m = 1: energy and period (small angle: 2 \pi / \sqrt(g))
  {
    pil_world *w = pil_world_new();
    pil_lsys *s = pil_add_lagrangian(w, 1, pend_T, pend_V, NULL);
    pil_lsys_set_state(s, (double[]){0.05}, (double[]){0});
    double T0 = 2 * M_PI / sqrt(9.81), hmax = 0, th_prev = 0.05, tcross = -1;
    double h0 = pil_lsys_jacobi(s);
    for (int i = 0; i < 4000; i++) {
      pil_world_step(w, 1e-3);
      hmax = fmax(hmax, fabs(pil_lsys_jacobi(s) - h0));
      double th = pil_lsys_q(s)[0];
      if (tcross < 0 && th_prev > 0 && th <= 0 && i > 10) {
        tcross = (i + 1) * 1e-3 + 0; // first zero crossing = T/4
      }
      th_prev = th;
    }
    check("pendulum energy drift", hmax, 0.0, 1e-8);
    check("pendulum quarter period", tcross, T0 / 4, 2e-3);
    pil_world_free(w);
  }
  // rotating arm: block sits at x = k L0 / (k - m Om^2) with no pendulum motion
  // (centrifugal term); stays put for ever
  {
    tower_params p = {1.0, 0.0001, 0.5, 60.0, 0.8, 3.0, 9.81};
    pil_world *w = pil_world_new();
    pil_lsys *s = pil_add_lagrangian(w, 2, tower_T, tower_V, &p);
    double xe = p.k * p.L0 / (p.k - (p.m1 + p.m2) * p.Om * p.Om);
    pil_lsys_set_state(s, (double[]){xe, 0}, (double[]){0, 0});
    for (int i = 0; i < 2000; i++) {
      pil_world_step(w, 1e-3);
    }
    check("centrifugal equilibrium", pil_lsys_q(s)[0], xe, 2e-4);
    pil_world_free(w);
  }
  // Full tower with Coriolis coupling: Jacobi integral conserved
  {
    tower_params p = {1.0, 0.6, 0.5, 60.0, 0.8, 2.0, 9.81};
    pil_world *w = pil_world_new();
    pil_lsys *s = pil_add_lagrangian(w, 2, tower_T, tower_V, &p);
    pil_lsys_set_state(s, (double[]){0.9, 0.7}, (double[]){0.3, -0.5});
    double h0 = pil_lsys_jacobi(s), dh = 0, xmin = 9, xmax = -9;
    for (int i = 0; i < 6000; i++) {
      pil_world_step(w, 1e-3);
      dh = fmax(dh, fabs(pil_lsys_jacobi(s) - h0));
      xmin = fmin(xmin, pil_lsys_q(s)[0]);
      xmax = fmax(xmax, pil_lsys_q(s)[0]);
    }
    check("tower Jacobi integral drift", dh, 0.0, 1e-5);
    check("tower block really moves", xmax - xmin > 0.05 ? 1.0 : 0.0, 1.0,
          1e-9);
    // Reset returns to initial state
    pil_world_reset(w);
    check("reset restores q", pil_lsys_q(s)[0], 0.9, 1e-12);
    pil_world_free(w);
  }
}

static void test_kirchhoff(void) {
  puts("Kirchhoff laws");
  pil_world *w = pil_world_new();
  pil_circuit *c = pil_add_circuit(w);
  pil_node g = pil_circuit_node(c, 0, 0), a = pil_circuit_node(c, 0, 2),
           x = pil_circuit_node(c, 2, 2), x2 = pil_circuit_node(c, 3, 2),
           y = pil_circuit_node(c, 3, 0), z = pil_circuit_node(c, 5, 2),
           zg = pil_circuit_node(c, 5, 0);
  pil_circuit_ground(c, g);
  pil_celement *vs = pil_vsource_dc(c, a, g, 12.0, "V");
  pil_celement *r1 = pil_resistor(c, a, x, 100.0, "R_1");
  pil_wire(c, x, x2);
  pil_celement *r2 = pil_resistor(c, x2, y, 300.0, "R_2");
  pil_celement *r3 = pil_resistor(c, x2, z, 150.0, "R_3");
  pil_wire(c, y, g);
  pil_wire(c, zg, g);
  pil_resistor(c, z, zg, 50.0, "R_4");
  pil_world_advance(w, 1e-3, 1e-4);

  // KCL at junction and at source's node
  check("KCL at junction", pil_circuit_node_current_sum(c, x), 0.0, 1e-9);
  check("KCL at source node", pil_circuit_node_current_sum(c, a), 0.0, 1e-9);
  // 2 parallel branches carry divider's currents: R3 + R4 = 200
  double rp = 1.0 / (1.0 / 300.0 + 1.0 / 200.0), itot = 12.0 / (100.0 + rp);
  check("series current", pil_celement_current(r1), itot, 1e-6);
  check("divider branch", pil_celement_current(r2), itot * rp / 300.0, 1e-6);
  check("parallel branch", pil_celement_current(r3), itot * rp / 200.0, 1e-6);
  // KVL around loop source - R1 - R2
  check("KVL loop",
        pil_celement_voltage(vs) -
            (pil_celement_voltage(r1) + pil_celement_voltage(r2)),
        0.0, 1e-6);
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
  test_link_and_plot();
  test_circuits();
  test_electromechanics();
  test_stops();
  test_circle_contact();
  test_lagrangian();
  test_kirchhoff();
  printf("\n%s (%d failure%s)\n", fails ? "FAILED" : "ALL PASSED", fails,
         fails == 1 ? "" : "s");
  return fails ? 1 : 0;
}
