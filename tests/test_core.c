/* Numerical regression tests against analytic solutions */

#include <math.h>
#include <pil/pil_biot.h>
#include <pil/pil_chaos.h>
#include <pil/pil_circuit.h>
#include <pil/pil_core.h>
#include <pil/pil_electrostatic.h>
#include <pil/pil_emwave.h>
#include <pil/pil_lagrange.h>
#include <pil/pil_nbody.h>
#include <pil/pil_ode.h>
#include <pil/pil_optics.h>
#include <pil/pil_orbit.h>
#include <pil/pil_particle.h>
#include <pil/pil_rotframe.h>
#include <pil/pil_scatter.h>
#include <pil/pil_solenoid.h>
#include <pil/pil_spherical.h>
#include <pil/pil_timer555.h>
#include <pil/pil_wave.h>
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

// 1. Mass on a spring: period 2*pi*sqrt(m/k), energy conserved
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

// 2. Rigid rod pendulum: small-angle period of a physical pendulum, energy
static void test_pendulum(void) {
  puts("physical pendulum (revolute)");
  pil_world *w = pil_world_new();
  double L = 1.0;
  double m = 3.0;
  double th0 = 0.05;
  // Uniform bar hinged at its end; COM at L/2
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

// 3. Cylinder rolling down an incline: a = g \sin(a) / 1.5
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

// 4. Rope: taut -> acts like a rod; goes slack when pushed together
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

// 6. Block resting on a floor (penalty contact): settles with N = m g
static void test_floor(void) {
  puts("block on floor (contact + friction)");
  pil_world *w = pil_world_new();
  pil_body *b = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.4, 0.2),
                                                 .mass = 5.0,
                                                 .pos = {0, 0.1},
                                                 .vel = {1.0, 0}});
  pil_add_floor_contact(w, b, pil_at(0, 0), (pil_vec2){1, 0}, 2e5, 2e3, 0.3);
  pil_world_advance(w, 0.2, 1e-4);
  // Sliding friction decelerates at \mu g = 2.943 m/s^2
  check("sliding deceleration", 1.0 - b->vel.x, 0.3 * 9.81 * 0.2, 2e-3);
  check("rest height (small penetration)", b->pos.y, 0.1, 2e-3);
  pil_world_free(w);
}

// 7. Double pendulum chain: constraint drift + energy over 10 s
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
// Rope goes slack if both are lifted
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

  // Equal masses at rest: tension = m g, nothing moves
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

// 9. Mouse grab: pick a body, drag its target; grabbed point follows
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

  // Grab at COM: pure translation, critically damped -> COM lands on target
  pil_world_set_grab(w, b, (pil_vec2){0, 0}, (pil_vec2){1.0, 0.5},
                     300.0 * b->mass, 2.0 * sqrt(300.0) * b->mass);
  pil_world_advance(w, 3.0, 1e-3);
  check("COM grab: x at target", b->pos.x, 1.0, 1e-3);
  check("COM grab: y at target", b->pos.y, 0.5, 1e-3);
  check("COM grab: no rotation", b->angle, 0.0, 1e-9);
  pil_world_clear_grab(w);

  // Off-centre grab: body rotates, but grabbed point reaches target
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
  // No slip: contact-point velocity along track is zero =>
  // \omega = -(D psidot)/r
  pil_world_advance(w, 0.3, 1e-3);
  pil_vec2 d = {ball->pos.x - c.x, ball->pos.y - c.y};
  double psid =
      (d.x * ball->vel.y - d.y * ball->vel.x) / (d.x * d.x + d.y * d.y);
  check("no slip: r*omega + D*psidot", r * ball->omega + D * psid, 0.0, 1e-6);
  pil_world_free(w);
}

// 11. Ball in a MOVING bowl: horizontal momentum and energy are conserved
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
                    (pil_vec2){1, 0}); // bowl slides horizontally
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

// 12. Ball over top of a cylinder (convex contact)
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

// 13. Weld: compound pendulum (bar + welded end mass), T from I_total/(m g d)
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

// 14. Rope snapping taut at speed: no visible stretch and snap is nelastic
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
  // Half-wave rectifier: peak output = V_pk - V_d with V_d from diode law
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
  // Switch closes at 1 ms, then RC charges
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
  // Motor spinning a flywheel against viscous drag:
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
  // Voice coil holding a spring-mounted mass: x_ss = K V / (R k)
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
  // Free mass hits a wall at x = 1 with restitution e:
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
  // Frictionless disc sliding inside a loop leaves wall where
  // \cos(\phi) = (2 g D - v0^2) / (3 g D), \phi measured from bottom
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
  // Rolling without slipping (high mu): bottom speed^2 = 4/3 g h for a solid
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
  double x = q[0];
  double xd = qd[0];
  double gd = qd[1];
  double sg = sin(q[1]);
  double cg = cos(q[1]);
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
    double T0 = 2 * M_PI / sqrt(9.81);
    double hmax = 0;
    double th_prev = 0.05;
    double tcross = -1;
    double h0 = pil_lsys_jacobi(s);
    for (int i = 0; i < 4000; i++) {
      pil_world_step(w, 1e-3);
      hmax = fmax(hmax, fabs(pil_lsys_jacobi(s) - h0));
      double th = pil_lsys_q(s)[0];
      if (tcross < 0 && th_prev > 0 && th <= 0 && i > 10) {
        tcross = (i + 1) * 1e-3 + 0; // First zero crossing = T/4
      }
      th_prev = th;
    }
    check("pendulum energy drift", hmax, 0.0, 1e-8);
    check("pendulum quarter period", tcross, T0 / 4, 2e-3);
    pil_world_free(w);
  }
  // Rotating arm: block sits at x = k L0 / (k - m Om^2) with no pendulum motion
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

typedef struct {
  double I1;
  double I3;
  double m;
  double l;
  double g;
} top_params;

static double top_T(const double *q, const double *qd, void *u) {
  const top_params *p = u;
  double st = sin(q[1]);
  double w3 = qd[2] + qd[0] * cos(q[1]);

  return 0.5 * p->I1 * (qd[1] * qd[1] + qd[0] * qd[0] * st * st) +
         0.5 * p->I3 * w3 * w3;
}

static double top_V(const double *q, void *u) {
  const top_params *p = u;

  return p->m * p->g * p->l * cos(q[1]);
}

static void test_top(void) {
  puts("heavy symmetric top (Euler angles)");
  top_params p = {0.03, 0.01, 1.0, 0.3, 9.81};
  pil_world *w = pil_world_new();
  pil_lsys *s = pil_add_lagrangian(w, 3, top_T, top_V, &p);
  pil_lsys_set_state(s, (double[]){2.3, 0.7, 0.0}, (double[]){4.2, 0.0, 70.0});
  const double *q = pil_lsys_q(s);
  const double *qd = pil_lsys_qd(s);
  double ppsi0 = p.I3 * (qd[2] + qd[0] * cos(q[1]));
  double pphi0 = p.I1 * qd[0] * sin(q[1]) * sin(q[1]) + ppsi0 * cos(q[1]);
  double h0 = pil_lsys_jacobi(s);
  double dpsi = 0;
  double dphi = 0;
  double dh = 0;
  for (int i = 0; i < 4000; i++) {
    pil_world_step(w, 1e-3);
    double ppsi = p.I3 * (qd[2] + qd[0] * cos(q[1]));
    double pphi = p.I1 * qd[0] * sin(q[1]) * sin(q[1]) + ppsi * cos(q[1]);
    dpsi = fmax(dpsi, fabs(ppsi - ppsi0));
    dphi = fmax(dphi, fabs(pphi - pphi0));
    dh = fmax(dh, fabs(pil_lsys_jacobi(s) - h0));
  }
  check("p_psi conserved", dpsi, 0.0, 1e-6);
  check("p_phi conserved", dphi, 0.0, 1e-6);
  check("energy conserved", dh, 0.0, 1e-5);
  pil_world_free(w);
}

static void test_rigid3(void) {
  puts("torque-free asymmetric body");
  pil_rigid3 p = {{1.0, 1.36, 1.64}};
  pil_world *w = pil_world_new();
  pil_lsys *s = pil_add_lagrangian(w, 3, pil_rigid3_T, NULL, &p);
  double w0[3] = {0.2, 3.0, 0.15};
  double q0[3], qd0[3];
  pil_rigid3_state(&p, w0, 0.0, q0, qd0);
  pil_lsys_set_state(s, q0, qd0);
  double Lm = sqrt(pow(p.I[0] * w0[0], 2) + pow(p.I[1] * w0[1], 2) +
                   pow(p.I[2] * w0[2], 2));
  double E0 = pil_lsys_jacobi(s);
  double dE = 0;
  double dL = 0;
  double dLz = 0;
  double dOrth = 0;
  double minw2 = 1e9;
  double maxw2 = -1e9;
  const double *q = pil_lsys_q(s);
  const double *qd = pil_lsys_qd(s);
  for (int i = 0; i < 20000; i++) {
    pil_world_step(w, 1e-3);
    double wb[3];
    pil_rigid3_omega(q, qd, wb);
    double Lb[3] = {p.I[0] * wb[0], p.I[1] * wb[1], p.I[2] * wb[2]};
    dE = fmax(dE, fabs(pil_lsys_jacobi(s) - E0));
    dL = fmax(dL,
              fabs(sqrt(Lb[0] * Lb[0] + Lb[1] * Lb[1] + Lb[2] * Lb[2]) - Lm));
    // L in space = R L_body must stay along z with length Lm
    double R[3][3];
    pil_rigid3_axes(q, R);
    double Lx = 0;
    double Ly = 0;
    double Lz = 0;
    for (int k = 0; k < 3; k++) {
      Lx += R[0][k] * Lb[k];
      Ly += R[1][k] * Lb[k];
      Lz += R[2][k] * Lb[k];
    }

    dLz = fmax(dLz, fmax(fabs(Lx), fmax(fabs(Ly), fabs(Lz - Lm))));
    double d = 0;
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        double dot = 0;
        for (int k = 0; k < 3; k++) {
          dot += R[k][a] * R[k][b];
        }

        d = fmax(d, fabs(dot - (a == b)));
      }
    }

    dOrth = fmax(dOrth, d);
    minw2 = fmin(minw2, wb[1]);
    maxw2 = fmax(maxw2, wb[1]);
  }

  check("energy conserved", dE, 0.0, 1e-5);
  check("|L| conserved", dL, 0.0, 1e-5);
  check("L stays fixed along space z", dLz, 0.0, 1e-5);
  check("axes stay orthonormal", dOrth, 0.0, 1e-12);
  // Spinning about intermediate axis 2 is unstable: omega_2 changes sign
  check("intermediate axis flips (w2 < 0 reached)", minw2 < 0 ? 1.0 : 0.0, 1.0,
        0);
  check("and returns (w2 > 2.5)", maxw2 > 2.5 ? 1.0 : 0.0, 1.0, 0);
  pil_world_free(w);
}

static void test_trace(void) {
  puts("traces");
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  pil_body *b = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.05),
                                                 .mass = 1.0,
                                                 .pos = {0, 0},
                                                 .vel = {2.0, 0}});
  pil_deco *tr = pil_add_trace(w, pil_on(b, 0, 0), 0.1, 100);
  pil_world_advance(w, 1.0, 1e-3);
  check("sample count", pil_trace_count(tr), 11, 1.5);
  check("first point", pil_trace_point(tr, 0).x, 0.0, 1e-9);
  check("last point", pil_trace_point(tr, pil_trace_count(tr) - 1).x, 2.0,
        0.03);
  pil_world_reset(w);
  check("cleared on reset", pil_trace_count(tr), 1, 0.5);
  pil_world_free(w);
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

static double plank_x(double t, void *u) {
  (void)u;

  return 0.5 * 2.0 * t * t; // a = 2 m/s^2
}

static void test_rolling_on_plank(void) {
  puts("disc rolling on a plank that is driven");
  // No slip on a plank accelerating at a: friction m a / 3, so a solid disc's
  // centre accelerates at a / 3 and contact point follows plank
  pil_world *w = pil_world_new();
  pil_body *pl = pil_add_body(
      w, &(pil_body_desc){.shape = pil_rect(2.0, 0.1), .mass = 4.0});
  pil_add_prismatic(w, pil_on(pl, 0, 0), pil_at(0, 0), (pil_vec2){1, 0});
  pil_add_coord_drive(w, pl, PIL_COORD_X, plank_x, NULL);
  pil_body *d = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(0.2),
                                                 .mass = 1.0,
                                                 .pos = {0.1, 0.25}});
  pil_add_rolling_on(w, d, (pil_vec2){0, 0}, 0.2, pil_on_at(pl, 0.1, 0.05),
                     (pil_vec2){1, 0});
  double x0 = d->pos.x;
  for (int i = 0; i < 1000; i++) {
    pil_world_step(w, 1e-3);
  }
  check("plank follows drive", pl->pos.x, plank_x(1.0, NULL), 1e-6);
  check("disc acceleration a/3", d->pos.x - x0, 0.5 * (2.0 / 3.0), 2e-4);
  check("disc stays on plank", d->pos.y, 0.25, 1e-6);
  check("constraint error", pil_world_constraint_error(w), 0.0, 1e-6);
  pil_world_free(w);
}

static void test_collision(void) {
  puts("disc collision, CM");
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  pil_body *a = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(0.3),
                                                 .mass = 1.0,
                                                 .pos = {-2, 0},
                                                 .vel = {2, 0}});
  pil_body *b = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(0.45),
                                                 .mass = 2.0,
                                                 .pos = {0.5, 0.2}});
  pil_add_collision(w, a, b, 0.3, 0.45, 0.8);
  pil_body *both[2] = {a, b};
  double T0 = pil_world_kinetic(w);
  pil_world_advance(w, 3.0, 1e-3);
  check("total px conserved", a->vel.x + 2.0 * b->vel.x, 2.0, 1e-9);
  check("total py conserved", a->vel.y + 2.0 * b->vel.y, 0.0, 1e-9);
  // Contact normal at impact: centres 0.75 apart, 0.2 offset across the path
  double c = sqrt(1.0 - (0.2 / 0.75) * (0.2 / 0.75));
  pil_vec2 n = {c, 0.2 / 0.75};
  double vn = (b->vel.x - a->vel.x) * n.x + (b->vel.y - a->vel.y) * n.y;
  // Approach speed along n is 2 cos(phi)
  check("separation speed = e * approach", vn, 0.8 * 2.0 * c, 1e-3);
  check("kinetic energy only decreases", pil_world_kinetic(w) < T0 ? 1.0 : 0.0,
        1.0, 0);
  pil_vec2 cm = pil_bodies_com(both, 2);
  check("CM moves uniformly (x)", cm.x,
        (-2.0 + 0.5 * 2.0) / 3.0 + 3.0 * 2.0 / 3.0, 1e-6);
  pil_world_free(w);

  // Equal masses, head on, e = 1: velocities are exchanged
  w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  a = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(0.2),
                                       .mass = 1.0,
                                       .pos = {-1, 0},
                                       .vel = {1, 0}});
  b = pil_add_body(
      w, &(pil_body_desc){.shape = pil_disc(0.2), .mass = 1.0, .pos = {1, 0}});
  pil_add_collision(w, a, b, 0.2, 0.2, 1.0);
  pil_world_advance(w, 3.0, 1e-3);
  check("elastic exchange: a at rest", a->vel.x, 0.0, 1e-9);
  check("elastic exchange: b takes v", b->vel.x, 1.0, 1e-9);
  pil_world_free(w);
}

static void test_central_force(void) {
  puts("central force: Kepler orbit");
  // GM = 1, perihelion 0.4, e = 0.6: a = 1, period 2 pi, areal velocity 0.4
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  pil_body *s = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(0.05),
                                                 .mass = 1.0,
                                                 .flags = PIL_STATIC});
  pil_body *p = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.03),
                                                 .mass = 0.01,
                                                 .pos = {0.4, 0},
                                                 .vel = {0, 2.0}});
  pil_add_gravity_pair(w, s, p, 1.0);
  double E0 = pil_world_energy(w);
  double rmax = 0;
  double amax = 0;
  double amin = 1e9;
  for (int i = 0; i < 6283; i++) {
    pil_world_advance(w, 1e-3, 1e-3);
    rmax = fmax(rmax, hypot(p->pos.x, p->pos.y));
    double rate = 0.5 * (p->pos.x * p->vel.y - p->pos.y * p->vel.x);
    amax = fmax(amax, rate);
    amin = fmin(amin, rate);
  }
  check("aphelion 2a - rp", rmax, 1.6, 2e-3);
  check("areal velocity constant", amax - amin, 0.0, 1e-6);
  check("energy conserved", (pil_world_energy(w) - E0) / fabs(E0), 0.0, 1e-6);
  check("back at perihelion x", p->pos.x, 0.4, 2e-3);
  pil_world_free(w);

  puts("Coulomb: like charges repel");
  w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  pil_body *a = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.03),
                                                 .mass = 1.0,
                                                 .pos = {0, 0}});
  pil_body *b = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.03),
                                                 .mass = 1.0,
                                                 .pos = {1, 0}});
  pil_add_coulomb(w, a, 1.0, b, 1.0, 2.0);
  pil_world_advance(w, 0.1, 1e-3);
  // F = k q q / r^2 = 2, a = 2, so after 0.1 s: dx = 0.5 * 2 * 0.01 = 0.01
  check("repel: b moves away", b->pos.x - 1.0, 0.01, 2e-4);
  check("repel: a moves back", a->pos.x, -0.01, 2e-4);
  pil_world_free(w);
}

static void test_lorentz_drag(void) {
  puts("charge in fields, drag");
  // q = m = B = 1, v = 1.5: circle of radius 1.5, period 2 pi, speed constant
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  pil_body *b = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.03),
                                                 .mass = 1.0,
                                                 .pos = {0, 0},
                                                 .vel = {1.5, 0}});
  pil_add_lorentz(w, b, 1.0, 0, 0, 1.0);
  pil_world_advance(w, 3.14159265358979, 1e-3);
  check("half turn: x back on axis", b->pos.x, 0.0, 2e-3);
  check("half turn: y = -2 r", b->pos.y, -3.0, 2e-3);
  check("speed constant", hypot(b->vel.x, b->vel.y), 1.5, 1e-4);
  pil_world_free(w);

  // Crossed fields, released at rest: mean drift E/B, cycloid height 2 E
  // m/(qB^2)
  w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  b = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.03),
                                       .mass = 1.0,
                                       .pos = {0, 0}});
  pil_add_lorentz(w, b, 1.0, 0, 0.5, 1.0);
  double E0 = pil_world_energy(w);
  double ymax = 0;
  for (int i = 0; i < 6283; i++) {
    pil_world_advance(w, 1e-3, 1e-3);
    if (b->pos.y > ymax) {
      ymax = b->pos.y;
    }
  }
  check("E x B drift over one turn", b->pos.x / 6.283185307, 0.5, 2e-3);
  check("cycloid height", ymax, 1.0, 2e-3);
  check("energy (kinetic + q E.x potential)", pil_world_energy(w), E0, 1e-4);
  pil_world_free(w);

  // Quadratic drag from rest-speed 10: v' = -c v^2 gives v = v0 / (1 + c v0 t)
  w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  b = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.03),
                                       .mass = 2.0,
                                       .vel = {10.0, 0}});
  pil_add_drag(w, b, 0.0, 0.4); // c = 0.2 per metre
  pil_world_advance(w, 1.0, 1e-3);
  check("quadratic drag speed", b->vel.x, 10.0 / (1.0 + 0.2 * 10.0), 1e-3);
  pil_world_free(w);
}

static void test_modes(void) {
  puts("normal modes");
  // k = 4, kc = 2, m = 1: omega1 = 2, omega2 = sqrt(8)
  // each mode returns after one of its own periods
  const double m = 1.0;
  const double k = 4.0;
  const double kc = 2.0;
  const double L = 1.2;
  double om[2] = {sqrt(k / m), sqrt((k + 2 * kc) / m)};
  double sgn[2] = {1.0, -1.0};
  for (int mode = 0; mode < 2; mode++) {
    pil_world *w = pil_world_new();
    pil_world_set_gravity(w, 0, 0);
    pil_body *b[2];
    for (int i = 0; i < 2; i++) {
      double x = (i + 1) * L + (i ? sgn[mode] : 1.0) * 0.1;
      b[i] = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.05),
                                              .mass = m,
                                              .pos = {x, 0}});
      pil_add_prismatic(w, pil_on(b[i], 0, 0), pil_at((i + 1) * L, 0),
                        (pil_vec2){1, 0});
    }
    pil_add_spring(w, pil_at(0, 0), pil_on(b[0], 0, 0), k, L);
    pil_add_spring(w, pil_on(b[0], 0, 0), pil_on(b[1], 0, 0), kc, L);
    pil_add_spring(w, pil_on(b[1], 0, 0), pil_at(3 * L, 0), k, L);
    double T = 2 * 3.14159265358979 / om[mode];
    pil_world_advance(w, T, 1e-4);
    check(mode ? "mode 2 returns after its period"
               : "mode 1 returns after its period",
          b[0]->pos.x, L + 0.1, 1e-4);
    pil_world_advance(w, 0.5 * T, 1e-4);
    check("and is at the opposite extreme after half more", b[0]->pos.x,
          L - 0.1, 1e-4);
    pil_world_free(w);
  }
}

static void test_bead_chain(void) {
  puts("bead chain standing waves");
  // N = 9 beads, zero-rest springs: mode n has omega = 2 sqrt(k/m) sin(n pi /
  // 2(N+1))
  const int N = 9;
  const double k = 100.0;
  const double m = 1.0;
  const double dx = 0.5;
  for (int mode = 1; mode <= 3; mode++) {
    pil_world *w = pil_world_new();
    pil_world_set_gravity(w, 0, 0);
    pil_body *b[9];
    for (int j = 0; j < N; j++) {
      double x = (j + 1) * dx;
      b[j] = pil_add_body(
          w, &(pil_body_desc){.shape = pil_point_mass(0.03),
                              .mass = m,
                              .pos = {x, 0.1 * sin(mode * 3.14159265358979 *
                                                   (j + 1) / (N + 1))}});
      pil_add_prismatic(w, pil_on(b[j], 0, 0), pil_at(x, 0), (pil_vec2){0, 1});
    }

    pil_add_spring(w, pil_at(0, 0), pil_on(b[0], 0, 0), k, 0.0);
    for (int j = 0; j + 1 < N; j++) {
      pil_add_spring(w, pil_on(b[j], 0, 0), pil_on(b[j + 1], 0, 0), k, 0.0);
    }

    pil_add_spring(w, pil_on(b[N - 1], 0, 0), pil_at((N + 1) * dx, 0), k, 0.0);
    double om =
        2.0 * sqrt(k / m) * sin(mode * 3.14159265358979 / (2.0 * (N + 1)));
    double T = 2.0 * 3.14159265358979 / om;
    pil_world_advance(w, T, 1e-4);
    double err = 0;
    for (int j = 0; j < N; j++) {
      err = fmax(err, fabs(b[j]->pos.y - 0.1 * sin(mode * 3.14159265358979 *
                                                   (j + 1) / (N + 1))));
    }

    check(mode == 1   ? "mode 1 returns after 2 pi / omega_1"
          : mode == 2 ? "mode 2 returns after 2 pi / omega_2"
                      : "mode 3 returns after 2 pi / omega_3",
          err, 0.0, 2e-4);
    pil_world_free(w);
  }
}

static void test_roll_race(void) {
  puts("rolling race: a = g sin(alpha) / (1 + I / m R^2)");
  const double al = 25.0 * 3.14159265358979 / 180.0;
  const double R = 0.2;
  const double g = 9.81;
  const double kI[3] = {1.0, 0.5, 0.4};
  for (int i = 0; i < 3; i++) {
    pil_world *w = pil_world_new();
    pil_vec2 dir = {cos(al), -sin(al)}, nrm = {sin(al), cos(al)};
    pil_body *b =
        pil_add_body(w, &(pil_body_desc){.shape = pil_disc(R),
                                         .mass = 1.0,
                                         .inertia = kI[i] * R * R,
                                         .pos = {nrm.x * R, nrm.y * R}});
    pil_add_rolling(w, b, (pil_vec2){0, 0}, R, (pil_vec2){0, 0}, dir);
    pil_world_advance(w, 1.0, 1e-3);
    double s = b->pos.x * dir.x + b->pos.y * dir.y -
               (nrm.x * dir.x + nrm.y * dir.y) * R;
    check(i == 0   ? "hoop distance in 1 s"
          : i == 1 ? "disc distance in 1 s"
                   : "sphere distance in 1 s",
          s, 0.5 * g * sin(al) / (1.0 + kI[i]), 1e-3);
    pil_world_free(w);
  }
}

static void res_push(pil_world *w, double t, void *user) {
  (void)w;
  struct {
    pil_body *b;
    double F0;
    double w;
  } *d = user;
  pil_apply_force(d->b, d->b->pos, pil_v2(d->F0 * cos(d->w * t), 0));
}

static void test_bridge_rectifier(void) {
  puts("bridge rectifier: v_load stays positive, ripple small, ~ peak - 2 "
       "diode drops");
  pil_world *w = pil_world_new();
  pil_circuit *c = pil_add_circuit(w);
  pil_node a = pil_circuit_node(c, 0, 1), b = pil_circuit_node(c, 0, 0),
           p = pil_circuit_node(c, 2, 2), n = pil_circuit_node(c, 2, -1);
  pil_circuit_ground(c, n);
  pil_vsource_ac(c, a, b, 12.0, 400.0, 0.0, "v");
  pil_diode(c, a, p, 1e-12, 0.04, "");
  pil_diode(c, b, p, 1e-12, 0.04, "");
  pil_diode(c, n, a, 1e-12, 0.04, "");
  pil_diode(c, n, b, 1e-12, 0.04, "");
  pil_celement *r = pil_resistor(c, p, n, 500.0, "R");
  pil_capacitor(c, p, n, 20e-6, 0.0, "C");
  pil_world_advance(w, 10e-3, 1e-4);
  double lo = 1e9;
  double hi = -1e9;
  for (int i = 0; i < 100; i++) {
    pil_world_advance(w, 5e-5, 1e-4);
    double v = pil_celement_voltage(r);
    if (v < lo) {
      lo = v;
    }
    if (v > hi) {
      hi = v;
    }
  }
  check("no negative load voltage", lo > 0.0, 1, 0);
  check("ripple under 1.5 V", hi - lo < 1.5, 1, 0);
  check("peak near 12 V - 2 diode drops (n_vt ln(I/Is) each)", hi,
        12.0 - 2 * 0.04 * log(0.03 / 1e-12), 1.0);
  pil_world_free(w);
}

static void test_resonance(void) {
  puts("driven oscillator at resonance: amplitude F0 / (c w0)");
  const double m = 1.0, k = 100.0, c = 2.0, F0 = 1.0, w0 = sqrt(k / m);
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  pil_body *b = pil_add_body(
      w,
      &(pil_body_desc){.shape = pil_rect(0.5, 0.4), .mass = m, .pos = {0, 0}});
  pil_add_prismatic(w, pil_on(b, 0, 0), pil_at(0, 0), (pil_vec2){1, 0});
  pil_add_spring(w, pil_at(-1, 0), pil_on(b, 0, 0), k, 1.0);
  pil_add_damper(w, pil_at(-1, 0), pil_on(b, 0, 0), c);
  struct {
    pil_body *b;
    double F0;
    double w;
  } d = {b, F0, w0};
  pil_add_custom_force(w, res_push, &d);
  pil_world_advance(w, 8.0, 5e-4);
  double amp = 0;
  for (int i = 0; i < 2000; i++) {
    pil_world_advance(w, 1e-3, 5e-4);
    if (fabs(b->pos.x) > amp)
      amp = fabs(b->pos.x);
  }
  check("resonant amplitude", amp, F0 / (c * w0), 5e-4);
  pil_world_free(w);
}

static void test_devices(void) {
  puts("transformer, op-amp, bjt, zener");
  // Transformer, turns ratio 2, loaded
  {
    pil_world *w = pil_world_new();
    pil_circuit *c = pil_add_circuit(w);
    pil_node p1 = pil_circuit_node(c, 0, 1), g = pil_circuit_node(c, 0, 0),
             s1 = pil_circuit_node(c, 3, 1);
    pil_circuit_ground(c, g);
    pil_vsource_ac(c, p1, g, 10.0, 1000.0, 0.0, "v");
    pil_celement *l1 = pil_inductor(c, p1, g, 1.0, 0.0, "L1");
    pil_celement *l2 = pil_inductor(c, s1, g, 4.0, 0.0, "L2");
    pil_celement *rl = pil_resistor(c, s1, g, 5e3, "R");
    pil_mutual(l1, l2, 0.999);
    pil_circuit_set_step(c, 1e-6);
    pil_world_advance(w, 5e-3, 1e-5);
    double pk = 0;
    for (int i = 0; i < 1000; i++) {
      pil_world_advance(w, 1e-6, 1e-5);
      double v = fabs(pil_celement_voltage(rl));
      if (v > pk)
        pk = v;
    }
    check("transformer secondary peak (n = 2)", pk, 20.0, 0.6);
    pil_world_free(w);
  }
  // Inverting op-amp: gain -Rf/Rin, saturates at the rail
  for (int k = 0; k < 2; k++) {
    pil_world *w = pil_world_new();
    pil_circuit *c = pil_add_circuit(w);
    pil_node in = pil_circuit_node(c, 0, 1), sum = pil_circuit_node(c, 1, 1),
             out = pil_circuit_node(c, 3, 1), g = pil_circuit_node(c, 0, 0);
    pil_circuit_ground(c, g);
    pil_vsource_dc(c, in, g, k ? 2.0 : 1.0, "v");
    pil_resistor(c, in, sum, 1e3, "Rin");
    pil_resistor(c, sum, out, 10e3, "Rf");
    pil_opamp(c, g, sum, out, 12.0, "A");
    pil_world_advance(w, 1e-3, 1e-4);
    check(k ? "op-amp saturates at -12 V" : "inverting gain -10",
          pil_circuit_node_voltage(c, out), k ? -12.0 : -10.0, 2e-3);
    pil_world_free(w);
  }
  // npn common emitter: active (Ic = beta Ib) and saturated
  for (int k = 0; k < 2; k++) {
    pil_world *w = pil_world_new();
    pil_circuit *c = pil_add_circuit(w);
    pil_node vcc = pil_circuit_node(c, 0, 3), col = pil_circuit_node(c, 2, 2),
             bas = pil_circuit_node(c, 1, 1), em = pil_circuit_node(c, 2, 0),
             vb = pil_circuit_node(c, 0, 1), g = pil_circuit_node(c, 0, 0);
    pil_circuit_ground(c, g);
    pil_vsource_dc(c, vcc, g, 10.0, "Vcc");
    pil_vsource_dc(c, vb, g, 5.0, "Vb");
    pil_resistor(c, vcc, col, 1e3, "Rc");
    pil_resistor(c, vb, bas, k ? 10e3 : 430e3, "Rb");
    pil_wire(c, em, g);
    pil_cdevice *q = pil_npn(c, col, bas, em, 100.0, "Q");
    pil_world_advance(w, 1e-3, 1e-4);
    if (!k) {
      check("npn Ic / Ib = beta",
            pil_cdevice_current(q, 0) / pil_cdevice_current(q, 1), 100.0, 5.0);
      check("npn active Vce", pil_cdevice_voltage(q, 0, 2), 9.0, 0.2);
    } else {
      check("npn saturated Vce < 0.3 V", pil_cdevice_voltage(q, 0, 2) < 0.3, 1,
            0);
    }

    check("KCL at base", pil_circuit_node_current_sum(c, bas), 0.0, 1e-9);
    pil_world_free(w);
  }
  // Zener shunt regulator
  for (int k = 0; k < 2; k++) {
    pil_world *w = pil_world_new();
    pil_circuit *c = pil_add_circuit(w);
    pil_node a = pil_circuit_node(c, 0, 2), z = pil_circuit_node(c, 2, 2),
             g = pil_circuit_node(c, 0, 0);
    pil_circuit_ground(c, g);
    pil_vsource_dc(c, a, g, k ? 15.0 : 12.0, "V");
    pil_resistor(c, a, z, 470.0, "Rs");
    pil_zener(c, g, z, 5.1, "Z");
    pil_resistor(c, z, g, 2e3, "RL");
    pil_world_advance(w, 1e-3, 1e-4);
    double vz = pil_circuit_node_voltage(c, z);
    check(k ? "zener holds at 15 V input" : "zener holds at 12 V input", vz,
          5.25, 0.2);
    pil_world_free(w);
  }
}

typedef struct {
  double m;
  double R;
  double Om;
  double g;
} hoop_t;
static double hoop_T_(const double *q, const double *qd, void *u) {
  const hoop_t *p = u;
  double s = sin(q[0]);
  return 0.5 * p->m * p->R * p->R * (qd[0] * qd[0] + p->Om * p->Om * s * s);
}
static double hoop_V_(const double *q, void *u) {
  const hoop_t *p = u;
  return -p->m * p->g * p->R * cos(q[0]);
}

static void test_hoop(void) {
  puts("bead on rotating hoop: equilibrium, small-oscillation frequency, "
       "Jacobi integral");
  hoop_t p = {1.0, 1.0, 5.0, 9.81};
  double th0 = acos(p.g / (p.R * p.Om * p.Om));
  pil_world *w = pil_world_new();
  pil_lsys *s = pil_add_lagrangian(w, 1, hoop_T_, hoop_V_, &p);
  pil_lsys_set_state(s, (double[]){th0 + 0.01}, (double[]){0});
  double h0 = pil_lsys_jacobi(s);
  // Timing first return to start: period 2 pi / (Om sin th0)
  double per = 2.0 * 3.14159265358979 / (p.Om * sin(th0));
  double tprev = 0;
  double tcross = 0;
  double qprev = pil_lsys_q(s)[0] - th0;
  int up = 0;
  for (int i = 0; i < 4000 && tcross == 0; i++) {
    pil_world_advance(w, 1e-3, 1e-4);
    double qn = pil_lsys_q(s)[0] - th0;
    if (qprev > 0 && qn <= 0) {
      if (up) {
        tcross = pil_world_time(w);
      } else {
        tprev = pil_world_time(w);
        up = 1;
      }
    }
    qprev = qn;
  }
  (void)tprev;
  check("equilibrium angle (g / R Om^2)", cos(th0), p.g / (p.R * p.Om * p.Om),
        1e-12);
  check("small-oscillation period", tcross - tprev, per, 2e-3);
  check("Jacobi integral conserved", pil_lsys_jacobi(s) - h0, 0.0, 1e-6);
  pil_world_free(w);
}

static double dp_T_(const double *q, const double *qd, void *u) {
  (void)u;
  return 0.5 * 2.0 * qd[0] * qd[0] + 0.5 * qd[1] * qd[1] +
         qd[0] * qd[1] * cos(q[0] - q[1]);
}
static double dp_V_(const double *q, void *u) {
  (void)u;
  return -2.0 * 9.81 * cos(q[0]) - 9.81 * cos(q[1]);
}

static void test_double_pendulum(void) {
  puts("double pendulum: energy conservation, small-angle normal modes");
  pil_world *w = pil_world_new();
  pil_lsys *s = pil_add_lagrangian(w, 2, dp_T_, dp_V_, NULL);
  pil_lsys_set_state(s, (double[]){2.0, 2.2}, (double[]){0, 0});
  double E0 = pil_lsys_jacobi(s);
  pil_world_advance(w, 10.0, 5e-4);
  check("energy drift over 10 s", pil_lsys_jacobi(s) - E0, 0.0, 1e-6);
  pil_world_free(w);
  // m1 = m2 = l1 = l2 = 1: slow normal mode has w^2 = (2 - sqrt 2) g and
  // theta_2 / theta_1 = sqrt 2; after one period the start state returns
  w = pil_world_new();
  s = pil_add_lagrangian(w, 2, dp_T_, dp_V_, NULL);
  double a = 1e-3;
  double r = sqrt(2.0);
  pil_lsys_set_state(s, (double[]){a, a * r}, (double[]){0, 0});
  double T = 2.0 * 3.14159265358979 / sqrt((2.0 - r) * 9.81);
  pil_world_advance(w, T, 2e-4);
  check("slow mode returns after one period (theta_1 / a)",
        pil_lsys_q(s)[0] / a, 1.0, 2e-3);
  check("slow mode keeps theta_2 / theta_1 = sqrt 2",
        pil_lsys_q(s)[1] / pil_lsys_q(s)[0], r, 2e-3);
  pil_world_free(w);
}

typedef struct {
  double m;
  double l;
  double g;
} sph_t;
static double sph_T(const double *q, const double *qd, void *u) {
  const sph_t *p = u;
  double s = sin(q[0]);
  return 0.5 * p->m * p->l * p->l * (qd[0] * qd[0] + s * s * qd[1] * qd[1]);
}
static double sph_V(const double *q, void *u) {
  const sph_t *p = u;
  return -p->m * p->g * p->l * cos(q[0]);
}

static void test_spherical(void) {
  puts("spherical pendulum: L_z and energy conserved, conical pendulum is "
       "steady");
  sph_t p = {1.0, 1.0, 9.81};
  pil_world *w = pil_world_new();
  pil_lsys *s = pil_add_lagrangian(w, 2, sph_T, sph_V, &p);
  pil_lsys_set_state(s, (double[]){0.9, 0}, (double[]){0, 2.2});
  double lz0 = sin(0.9) * sin(0.9) * 2.2, e0 = pil_lsys_jacobi(s);
  pil_world_advance(w, 20.0, 5e-4);
  double st = sin(pil_lsys_q(s)[0]);
  check("L_z drift", st * st * pil_lsys_qd(s)[1] - lz0, 0.0, 1e-6);
  check("energy drift", pil_lsys_jacobi(s) - e0, 0.0, 1e-6);
  pil_world_free(w);
  // Circular motion: phi'^2 = g / (l cos th)
  w = pil_world_new();
  s = pil_add_lagrangian(w, 2, sph_T, sph_V, &p);
  double th = 0.7;
  pil_lsys_set_state(s, (double[]){th, 0}, (double[]){0, sqrt(9.81 / cos(th))});
  pil_world_advance(w, 10.0, 5e-4);
  check("conical pendulum keeps theta", pil_lsys_q(s)[0], th, 1e-6);
  pil_world_free(w);
}

static int doppler_count(double xo, double t) {
  int n = 0;
  for (int i = 0; i * 0.6 <= t; i++) {
    double te = i * 0.6, xs = -3.0 + 0.5 * te;
    if (te + fabs(xo - xs) <= t) {
      n++;
    }
  }
  return n;
}

static void test_doppler(void) {
  puts("doppler: reception rates f / (1 -+ v/c)");
  double f = 1.0 / 0.6;
  double ra = (doppler_count(5.5, 16.0) - doppler_count(5.5, 10.0)) /
              6.0; // Source passes A at t = 17
  double rb = (doppler_count(-5.5, 60.0) - doppler_count(-5.5, 30.0)) / 30.0;
  check("observer ahead: f / (1 - v/c)", ra, f / 0.5, 0.4);
  check("observer behind: f / (1 + v/c)", rb, f / 1.5, 0.05);
}

static void test_mosfet(void) {
  puts("mosfet: square law, cmos inverter");
  // Common-Source NMOS: Id = K/2 (Vgs - Vth)^2 (1 + 0.02 (Vds - Vov))
  {
    pil_world *w = pil_world_new();
    pil_circuit *c = pil_add_circuit(w);
    pil_node vdd = pil_circuit_node(c, 0, 3), d = pil_circuit_node(c, 2, 2),
             g = pil_circuit_node(c, 0, 1), z = pil_circuit_node(c, 0, 0);
    pil_circuit_ground(c, z);
    pil_vsource_dc(c, vdd, z, 5.0, "Vdd");
    pil_vsource_dc(c, g, z, 2.0, "Vg");
    pil_resistor(c, vdd, d, 1e3, "R");
    pil_cdevice *m = pil_nmos(c, d, g, z, 2e-3, 1.0, "M");
    pil_world_advance(w, 1e-3, 1e-4);
    double vds = pil_circuit_node_voltage(c, d), id = pil_cdevice_current(m, 0);
    check("saturation current", id, 1e-3 * (1.0 + 0.02 * (vds - 1.0)), 1e-9);
    check("KCL: drain resistor current = Id", (5.0 - vds) / 1e3, id, 1e-9);
    check("gate draws no current", pil_cdevice_current(m, 1), 0.0, 1e-15);
    pil_world_free(w);
  }
  // CMOS inverter: low in -> high out, high in -> low out, midpoint at Vdd / 2
  double vin[3] = {0.0, 5.0, 2.5}, want[3] = {5.0, 0.0, 2.5};
  for (int k = 0; k < 3; k++) {
    pil_world *w = pil_world_new();
    pil_circuit *c = pil_add_circuit(w);
    pil_node vdd = pil_circuit_node(c, 2, 4), out = pil_circuit_node(c, 2, 2),
             in = pil_circuit_node(c, 0, 2), z = pil_circuit_node(c, 0, 0);
    pil_circuit_ground(c, z);
    pil_vsource_dc(c, vdd, z, 5.0, "Vdd");
    pil_vsource_dc(c, in, z, vin[k], "Vin");
    pil_pmos(c, out, in, vdd, 2e-3, 1.0, "Mp");
    pil_nmos(c, out, in, z, 2e-3, 1.0, "Mn");
    pil_resistor(c, out, z, 1e9, "Rl");
    pil_world_advance(w, 1e-3, 1e-4);
    check(k == 0   ? "inverter in=0 -> 5 V"
          : k == 1 ? "inverter in=5 -> 0 V"
                   : "inverter in=2.5 -> 2.5 V",
          pil_circuit_node_voltage(c, out), want[k], 0.05);
    pil_world_free(w);
  }
}

static double sine_f0(double t, void *u) {
  (void)u;
  return sin(2.0 * 3.14159265358979 * 503.29212 * t);
}

static void test_rlc_response(void) {
  puts("series RLC at resonance: gain 1 across R, Q across C");
  pil_world *w = pil_world_new();
  pil_circuit *c = pil_add_circuit(w);
  pil_node a = pil_circuit_node(c, 0, 1), n1 = pil_circuit_node(c, 1, 1),
           n2 = pil_circuit_node(c, 2, 1), g = pil_circuit_node(c, 0, 0);
  pil_circuit_ground(c, g);
  pil_vsource_fn(c, a, g, sine_f0, NULL, "v");
  pil_celement *r = pil_resistor(c, a, n1, 10.0, "R");
  pil_inductor(c, n1, n2, 10e-3, 0.0, "L");
  pil_celement *cap = pil_capacitor(c, n2, g, 10e-6, 0.0, "C");
  pil_circuit_set_step(c, 2e-6);
  pil_world_advance(w, 30e-3, 1e-4);
  double vr = 0, vc = 0;
  for (int i = 0; i < 4000; i++) {
    pil_world_advance(w, 1e-6 * 5, 1e-4);
    double a1 = fabs(10.0 * pil_celement_current(r)),
           a2 = fabs(pil_celement_voltage(cap));
    if (a1 > vr)
      vr = a1;
    if (a2 > vc)
      vc = a2;
  }
  check("band-pass gain at f0", vr, 1.0, 0.01);
  check("low-pass gain at f0 = Q", vc, sqrt(10e-3 / 10e-6) / 10.0, 0.03);
  pil_world_free(w);
}

static void uniform_b(const double x[3], double t, double B[3], double E[3],
                      void *u) {
  (void)x;
  (void)t;
  (void)u;
  B[0] = B[1] = 0;
  B[2] = 2.0;
  E[0] = E[1] = E[2] = 0;
}

static void bottle_b(const double x[3], double t, double B[3], double E[3],
                     void *u) {
  (void)t;
  (void)u;
  B[0] = -3.0 * x[0] * x[2] / 9.0;
  B[1] = -3.0 * x[1] * x[2] / 9.0;
  B[2] = 3.0 * (1.0 + x[2] * x[2] / 9.0);
  E[0] = E[1] = E[2] = 0;
}

static void test_boris(void) {
  puts("charged particle pusher: gyration and magnetic bottle");
  pil_particle3 p = {{1, 0, 0}, {0, 2.0, 0.5}, 1.0};
  double v0 = sqrt(4.0 + 0.25);
  // omega = qB/m = 2: after one period 2 pi / 2 the particle is back in x, y
  int n = 3142;
  double dt = 3.14159265358979 / n;
  for (int i = 0; i < n; i++)
    pil_particle3_step(&p, uniform_b, NULL, 0, dt);
  check("speed conserved",
        sqrt(p.v[0] * p.v[0] + p.v[1] * p.v[1] + p.v[2] * p.v[2]), v0, 1e-12);
  check("gyro-radius circle closes (x)", p.x[0], 1.0, 2e-4);
  check("drift along B (z = v_par T)", p.x[2], 0.5 * 3.14159265358979, 1e-6);
  // Bottle: mu = v_perp^2 / (2 B) conserved, mirror point z_m = L sqrt(1/sin^2
  // - 1)
  pil_particle3 q = {{0.33, 0, 0}, {0, 1.0, 0.35}, 1.0};
  double mu0 = 0;
  double mumax = 0;
  double mumin = 1e9;
  double zmax = 0;
  for (int i = 0; i < 15000; i++) {
    double B[3];
    double E[3];
    bottle_b(q.x, 0, B, E, NULL);
    double bm = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
    double vpar = (q.v[0] * B[0] + q.v[1] * B[1] + q.v[2] * B[2]) / bm;
    double v2 = q.v[0] * q.v[0] + q.v[1] * q.v[1] + q.v[2] * q.v[2];
    double mu = 0.5 * (v2 - vpar * vpar) / bm;
    if (i == 0) {
      mu0 = mu;
    }
    if (mu > mumax) {
      mumax = mu;
    }
    if (mu < mumin) {
      mumin = mu;
    }
    if (fabs(q.x[2]) > zmax) {
      zmax = fabs(q.x[2]);
    }
    pil_particle3_step(&q, bottle_b, NULL, 0, 0.004);
  }
  check("mu varies by less than 2 %", (mumax - mumin) / mu0 < 0.02, 1, 0);
  check("mirror point z_m = L sqrt(v_par^2 / v_perp^2)", zmax, 3.0 * 0.35,
        0.08);
}

static void dipole_b(const double x[3], double t, double B[3], double E[3],
                     void *u) {
  (void)t;
  (void)u;
  double r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  double r5 = r2 * r2 * sqrt(r2);
  B[0] = -3.0 * 100.0 * x[2] * x[0] / r5;
  B[1] = -3.0 * 100.0 * x[2] * x[1] / r5;
  B[2] = 100.0 * (r2 - 3.0 * x[2] * x[2]) / r5;
  E[0] = E[1] = E[2] = 0;
}

static void test_dipole(void) {
  puts("dipole field: trapped proton bounces, drifts west, mu conserved");
  pil_particle3 p = {{2.79, 0, 0}, {0, 1.0, 0.5}, 1.0};
  double v0 = sqrt(1.25);
  double zmax = 0;
  double mumin = 1e9;
  double mumax = 0;
  double mu0 = 0;
  double phi = 0;
  double prev = 0;
  for (int i = 0; i < 25000; i++) {
    double B[3];
    double E[3];
    dipole_b(p.x, 0, B, E, NULL);
    double bm = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
    double vpar = (p.v[0] * B[0] + p.v[1] * B[1] + p.v[2] * B[2]) / bm;
    double v2 = p.v[0] * p.v[0] + p.v[1] * p.v[1] + p.v[2] * p.v[2];
    double mu = 0.5 * (v2 - vpar * vpar) / bm;
    if (i == 0) {
      mu0 = mu;
    }
    if (i > 0 && i % 25 == 0) { // Average over the gyration
      if (mu > mumax) {
        mumax = mu;
      }
      if (mu < mumin) {
        mumin = mu;
      }
    }
    if (fabs(p.x[2]) > zmax) {
      zmax = fabs(p.x[2]);
    }
    double ph = atan2(p.x[1], p.x[0]);
    double d = ph - prev;
    if (d > 3.14159265358979) {
      d -= 2 * 3.14159265358979;
    }
    if (d < -3.14159265358979) {
      d += 2 * 3.14159265358979;
    }
    if (i > 0) {
      phi += d;
    }
    prev = ph;
    pil_particle3_step(&p, dipole_b, NULL, 0, 0.004);
  }
  check("speed conserved",
        sqrt(p.v[0] * p.v[0] + p.v[1] * p.v[1] + p.v[2] * p.v[2]), v0, 1e-9);
  check("bounces: |z| stays within mirror points", zmax > 0.45 && zmax < 0.75,
        1, 0);
  check("proton drifts westward (longitude decreases)", phi < -2.0, 1, 0);
  (void)mu0;
  (void)mumin;
  (void)mumax;
}

static void test_foucault(void) {
  puts("Foucault pendulum: swing plane turns at -Omega sin(lat)");
  double W0 = 6.283185307179586, WR = 0.15, h = 0.002;
  double s[4] = {1, 0, 0, 0};
  double xs[30000];
  double ys[30000];
  for (int i = 0; i < 30000; i++) {
    xs[i] = s[0];
    ys[i] = s[1];
    double k[4][4];
    double u[4];
    for (int st = 0; st < 4; st++) {
      double c = st == 0 ? 0 : (st == 3 ? 1 : 0.5);
      for (int j = 0; j < 4; j++) {
        u[j] = s[j] + (st ? c * h * k[st - 1][j] : 0);
      }
      k[st][0] = u[2];
      k[st][1] = u[3];
      k[st][2] = -W0 * W0 * u[0] + 2 * WR * u[3];
      k[st][3] = -W0 * W0 * u[1] - 2 * WR * u[2];
    }
    for (int j = 0; j < 4; j++)
      s[j] += h / 6 * (k[0][j] + 2 * k[1][j] + 2 * k[2][j] + k[3][j]);
  }
  double sxx = 0;
  double syy = 0;
  double sxy = 0;
  for (int j = 29500 - 500; j < 29500; j++) { // One period centred at t = 58.5
    sxx += xs[j] * xs[j];
    syy += ys[j] * ys[j];
    sxy += xs[j] * ys[j];
  }
  double psi = 0.5 * atan2(2 * sxy, sxx - syy);
  double want = -WR * 58.5;
  want = fmod(want + 100 * 3.14159265358979,
              3.14159265358979); // Plane angle is mod pi
  psi = fmod(psi + 100 * 3.14159265358979, 3.14159265358979);
  check("plane angle at t=58.5 s equals -W t (mod pi)", psi, want, 0.03);
}

static void test_chaos(void) {
  puts("driven pendulum: sensitive dependence on initial conditions");
  double a[2] = {0.2, 0};
  double b[2] = {0.2 + 1e-9, 0};
  double T = 2 * 3.14159265358979 / (2.0 / 3.0), h = T / 600;
  double d0 = 1e-9;
  double d = 0;
  double tm = 0;
  for (int i = 0; i < 600 * 40; i++) {
    for (int k = 0; k < 2; k++) {
      double *s = k ? b : a;
      double k1[2];
      double k2[2];
      double k3[2];
      double k4[2];
      double u[2];
      double tt = tm;
#define F(o, in, t)                                                            \
  (o)[0] = (in)[1];                                                            \
  (o)[1] = -(in)[1] / 2.0 - sin((in)[0]) + 1.5 * cos((2.0 / 3.0) * (t));
      F(k1, s, tt)
      u[0] = s[0] + 0.5 * h * k1[0];
      u[1] = s[1] + 0.5 * h * k1[1];
      F(k2, u, tt + 0.5 * h)
      u[0] = s[0] + 0.5 * h * k2[0];
      u[1] = s[1] + 0.5 * h * k2[1];
      F(k3, u, tt + 0.5 * h)
      u[0] = s[0] + h * k3[0];
      u[1] = s[1] + h * k3[1];
      F(k4, u, tt + h)
#undef F
      for (int j = 0; j < 2; j++) {
        s[j] += h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
      }
    }
    tm += h;
  }
  d = fabs(a[0] - b[0]) + fabs(a[1] - b[1]);
  check("separation grew by > 1e4 in 40 drive periods", d > 1e4 * d0, 1, 0);
}

static void biot_pair_cb(double x, double z, double *bx, double *bz, void *u) {
  pil_biot_pair(*(double *)u, x, z, 400, bx, bz);
}

static void test_helmholtz(void) {
  puts("Biot-Savart: loop on-axis field, Helmholtz flatness, field lines");
  double bx;
  double bz;
  for (int k = 0; k < 3; k++) {
    double z = 0.4 * k - 0.4;
    pil_biot_loop(1.0, 0.1, 0.0, z, 400, &bx, &bz);
    check("single loop vs 2 pi r^2 / (r^2 + z^2)^1.5", bz,
          pil_biot_loop_axis(1.0, 0.1, z), 1e-6);
    check("no radial field on the axis", bx, 0.0, 1e-9);
  }

  double b0;
  double b1;
  double b2;
  double d = 1.0;
  pil_biot_pair(d, 0, 0, 400, &bx, &b0);
  pil_biot_pair(d, 0, 0.2, 400, &bx, &b1);
  pil_biot_pair(d, 0, -0.2, 400, &bx, &b2);
  check("Helmholtz: B(+-0.2) matches B(0) to < 0.5 % (flat)",
        fabs(b1 / b0 - 1) < 0.005 && fabs(b2 / b1 - 1) < 1e-9, 1, 0);
  d = 1.5;
  pil_biot_pair(d, 0, 0, 400, &bx, &b0);
  pil_biot_pair(d, 0, 0.2, 400, &bx, &b1);
  check("d = 1.5 R: field dips in the middle", b1 > b0, 1, 0);
  // A field line through axis stays on axis, directed along +z
  static double px[300];
  static double pz[300];
  int n =
      pil_field_line(biot_pair_cb, &d, 0.0, 0.0, 1, 0.05, 200, 5, 5, px, pz);
  check("axis field line is straight (x = 0)", fabs(px[n - 1]), 0.0, 1e-9);
  check("axis field line runs along +z", pz[n - 1] > 5.0 - 0.06, 1, 0);
}

// static void test_helmholtz(void) {
//   puts("Biot-Savart: loop on-axis field and Helmholtz flatness");
//   // Unit current, mu0/4pi = 1: loop of radius 1 gives B_z = 2 pi / (1 +
//   // z^2)^(3/2)
//   double d = 1.0;
//   double B[3];
//   for (int k = 0; k < 3; k++) {
//     double z = 0.2 * k - 0.2;
//     double tot = 0;
//     for (int l = 0; l < 2; l++) {
//       double z0 = l ? d / 2 : -d / 2;
//       double bz = 0;
//       for (int i = 0; i < 400; i++) {
//         double ph = 2 * 3.14159265358979 * (i + 0.5) / 400;
//         double dph = 2 * 3.14159265358979 / 400;
//         double px = cos(ph);
//         double py = sin(ph);
//         double dlx = -py * dph;
//         double dly = px * dph;
//         double rx = -px;
//         double ry = -py;
//         double rz = z - z0;
//         double r3 = pow(rx * rx + ry * ry + rz * rz, 1.5);
//         bz += (dlx * ry - dly * rx) / r3;
//       }
//       tot += bz;
//       check("single loop vs 2 pi / (1+z^2)^1.5", bz,
//             2 * 3.14159265358979 / pow(1 + (z - z0) * (z - z0), 1.5), 1e-6);
//     }
//     B[k] = tot;
//   }
//   check("Helmholtz: B(+-0.2) matches B(0) to < 0.5 % (flat)",
//         fabs(B[0] - 0.5 * (B[0] + B[2])) < 1e-9 &&
//             fabs(B[0] / B[1] - 1) < 0.005,
//         1, 0);
// }

static void test_hohmann(void) {
  puts("Hohmann transfer r1 = 1 -> r2 = 2: two burns, flight time, final "
       "circle");
  double a = 1.5;
  double tof = 3.14159265358979 * sqrt(a * a * a);
  double dt = 5e-4;
  double dv1 = (sqrt(2.0 * 2 / 3) - 1);
  double dv2 = sqrt(0.5) * (1 - sqrt(2.0 / 3));
  double x = 1;
  double y = 0;
  double vx = 0;
  double vy = 1;
  double rmax = 0;
  double rmin = 1e9;
  int b1 = 0;
  double b2 = (int)(tof / dt + 0.5);
  double n = b2 + 12000;
  for (int i = 0; i < n; i++) {
    double r = hypot(x, y);
    if (i == b1 || i == b2) {
      double sp = hypot(vx, vy);
      double dv = i == b1 ? dv1 : dv2;
      vx += dv * vx / sp;
      vy += dv * vy / sp;
    }
    if (i > b2 + 100 && r > rmax) {
      rmax = r;
    }
    if (i > b2 + 100 && r < rmin) {
      rmin = r;
    }
    double r3 = r * r * r;
    vx -= 0.5 * dt * x / r3;
    vy -= 0.5 * dt * y / r3;
    x += dt * vx;
    y += dt * vy;
    r = hypot(x, y);
    r3 = r * r * r;
    vx -= 0.5 * dt * x / r3;
    vy -= 0.5 * dt * y / r3;
  }
  check("final orbit is circular at r = 2 (rmax)", rmax, 2.0, 2e-3);
  check("final orbit is circular at r = 2 (rmin)", rmin, 2.0, 2e-3);
  check("specific energy = -1/(2 r2)",
        0.5 * (vx * vx + vy * vy) - 1.0 / hypot(x, y), -0.25, 1e-3);
}

static void test_nbody(void) {
  puts("n-body: figure-eight choreography, conservation laws, Kepler pair");
  pil_nbody s;
  pil_nbody_figure_eight(&s);
  double x0[3][2];
  double e0 = pil_nbody_energy(&s);
  double l0 = pil_nbody_angmom(&s);
  double p[2];
  double c[2];
  for (int b = 0; b < 3; b++) {
    x0[b][0] = s.x[b][0];
    x0[b][1] = s.x[b][1];
  }
  pil_nbody_momentum(&s, p);
  check("figure-eight has zero momentum", hypot(p[0], p[1]), 0.0, 1e-8);
  pil_nbody_com(&s, c);
  check("figure-eight centre of mass at origin", hypot(c[0], c[1]), 0.0, 1e-8);
  int n = (int)(PIL_FIGURE_EIGHT_PERIOD / 1e-4 + 0.5);
  double emax = 0;
  for (int i = 0; i < n; i++) {
    pil_nbody_step(&s, PIL_FIGURE_EIGHT_PERIOD / n);
    double de = fabs(pil_nbody_energy(&s) - e0);
    if (de > emax) {
      emax = de;
    }
  }

  check("energy error stays below 1e-6", emax < 1e-6, 1, 0);
  check("angular momentum conserved", pil_nbody_angmom(&s) - l0, 0.0, 1e-8);
  double back = 0;
  for (int b = 0; b < 3; b++) {
    back += hypot(s.x[b][0] - x0[b][0], s.x[b][1] - x0[b][1]);
  }
  check("every body is back after one period", back, 0.0, 2e-4);
  // Two bodies: equal masses on a circular mutual orbit, m = 1, separation 2,
  // v = 1/sqrt(2*2)... v^2 = G m / (2 d)
  pil_nbody k = {
      2, 1.0, 0.0, {1.0, 1.0}, {{-1, 0}, {1, 0}}, {{0, -0.5}, {0, 0.5}}};
  double ek = pil_nbody_energy(&k);
  for (int i = 0; i < 20000; i++) {
    pil_nbody_step(&k, 1e-3);
  }
  check("circular pair keeps separation 2",
        hypot(k.x[1][0] - k.x[0][0], k.x[1][1] - k.x[0][1]), 2.0, 1e-5);
  check("circular pair keeps energy", pil_nbody_energy(&k), ek, 1e-8);
}

static void test_emwave(void) {
  puts("EM plane wave: E x B along k, B = E / c, u_E = u_B, S = c u, radiation "
       "pressure");
  pil_emwave w = {2.0, 1.5, 2.0 * M_PI / 4.0, 0.0, 0.3};
  double E[3];
  double B[3];
  double S[3];
  double sum_s = 0;
  double smin = 1e9;
  int n = 4000;
  for (int i = 0; i < n; i++) {
    double t = 2.0 * i / n * 1.0; // One period T = lambda / c = 2
    double x = 0.37;
    pil_emwave_field(&w, x, t, E, B);
    pil_emwave_poynting(&w, x, t, S);
    double u = pil_emwave_energy_density(&w, x, t);
    if (i % 400 == 0) {
      check("B = E / c (magnitude)", hypot(B[1], B[2]), hypot(E[1], E[2]) / 2.0,
            1e-12);
      check("E . B = 0", E[1] * B[1] + E[2] * B[2], 0.0, 1e-12);
      check("S = c u along +x", S[0], 2.0 * u, 1e-12);
      check("S has no transverse part", hypot(S[1], S[2]), 0.0, 1e-12);
    }

    sum_s += S[0];
    if (S[0] < smin) {
      smin = S[0];
    }
  }

  check("S never negative", smin >= -1e-12, 1, 0);
  check("mean S = c e0^2 / 2", sum_s / n, pil_emwave_intensity(&w), 1e-9);
  check("intensity = c e0^2 / 2", pil_emwave_intensity(&w),
        0.5 * 2.0 * 1.5 * 1.5, 1e-12);
  check("reflector feels twice the pressure", pil_emwave_pressure(&w, 1),
        2.0 * pil_emwave_pressure(&w, 0), 1e-12);
  // Circular polarisation: constant energy density and flux
  pil_emwave c = {1.0, 1.0, 1.0, 1.0, 0.0};
  double u0 = pil_emwave_energy_density(&c, 0, 0);
  double u1 = pil_emwave_energy_density(&c, 0.7, 1.9);
  check("circular polarisation: u is constant", u1, u0, 1e-12);
  check("circular intensity = c e0^2", pil_emwave_intensity(&c), 1.0, 1e-12);
}

static void test_rutherford(void) {
  puts("Rutherford scattering: deflection, closest approach, cross section");
  static double path[20001][2];
  double K = 2.0, v0 = 2.0;
  for (int k = 0; k < 4; k++) {
    double b = 0.2 + 0.5 * k;
    double th;
    double rm;
    pil_scatter_trace(K, v0, b, -12.0, 0.002, 20000, path, &th, &rm);
    check("measured deflection = 2 atan(kappa / b)", th,
          pil_scatter_angle(K, v0, b), 1e-2);
    check("closest approach = kappa + sqrt(kappa^2 + b^2)", rm,
          pil_scatter_rmin(K, v0, b), 2e-3);
  }

  // Lower side: negative impact parameter deflects other way
  double th;
  pil_scatter_trace(K, v0, -0.7, -12.0, 0.002, 20000, path, &th, NULL);
  check("b < 0 deflects towards -y", th, -pil_scatter_angle(K, v0, 0.7), 1e-2);
  check("head-on limit: r_min = 2 kappa", pil_scatter_rmin(K, v0, 0.0),
        2.0 * K / (v0 * v0), 1e-12);
  check("b -> infinity: no deflection", pil_scatter_angle(K, v0, 1e6), 0.0,
        1e-5);
  // Impact parameter and cross section: dsigma = (b / sin th) |db / dth|
  double t0 = 1.1;
  double d = 1e-6;
  double bb = pil_scatter_impact(K, v0, t0),
         dbdt = (pil_scatter_impact(K, v0, t0 + d) -
                 pil_scatter_impact(K, v0, t0 - d)) /
                (2 * d);
  check("impact(angle(b)) = b", pil_scatter_angle(K, v0, bb), t0, 1e-12);
  check("dsigma = b |db/dtheta| / sin(theta)", pil_scatter_dsigma(K, v0, t0),
        bb * fabs(dbdt) / sin(t0), 1e-6);
}

static void test_electrostatic(void) {
  puts("electrostatics: Coulomb field, Gauss's law, field lines, energy");
  pil_charges c = {2, {1.0, -0.4}, {{-1.5, 0}, {1.5, 0}}};
  double ex;
  double ey;
  pil_charges one = {1, {2.0}, {{0, 0}}};
  pil_es_field(&one, 3.0, 4.0, &ex, &ey);
  check("point charge: |E| = q / r^2", hypot(ex, ey), 2.0 / 25.0, 1e-12);
  check("point charge: V = q / r", pil_es_potential(&one, 3.0, 4.0), 2.0 / 5.0,
        1e-12);
  // Gauss: flux through any sphere around both charges = 4 pi (q1 + q2); only
  // enclosed charge counts
  check("flux through R = 3 sphere = 4 pi (q1 + q2)",
        pil_es_flux_sphere(&c, 0.0, 3.0, 400), 4.0 * M_PI * 0.6, 1e-4);
  check("flux does not depend on the radius",
        pil_es_flux_sphere(&c, 0.3, 7.0, 400), 4.0 * M_PI * 0.6, 1e-4);
  check("flux through a sphere enclosing only q1 = 4 pi q1",
        pil_es_flux_sphere(&c, -1.5, 0.8, 400), 4.0 * M_PI * 1.0, 1e-4);
  check("charges outside give zero flux",
        pil_es_flux_sphere(&c, -9.0, 0.5, 400), 0.0, 1e-5);
  // Dipole: V = 0 on mid-plane, field there parallel to axis
  pil_charges d = {2, {1.0, -1.0}, {{-1, 0}, {1, 0}}};
  check("dipole: V = 0 on the mid-plane", pil_es_potential(&d, 0.0, 1.3), 0.0,
        1e-12);
  pil_es_field(&d, 0.0, 1.3, &ex, &ey);
  check("dipole: E on the mid-plane has no y part", ey, 0.0, 1e-12);
  check("dipole energy = -q^2 / r12", pil_es_energy(&d), -0.5, 1e-12);
  // A field line from + ends on - (opposite equal charges)
  static double px[4000], py[4000];
  int n = pil_es_field_line(&d, -1.0 + 0.12 * cos(2.0), 0.12 * sin(2.0), 1, 0,
                            0.02, 3999, 20, 20, 0.1, px, py);
  check("line from + ends at the - charge",
        hypot(px[n - 1] - 1.0, py[n - 1]) < 0.12, 1, 0);
  // Like charges: no line connects them, lines run off to infinity
  pil_charges l = {2, {1.0, 1.0}, {{-1, 0}, {1, 0}}};
  n = pil_es_field_line(&l, -1.0 + 0.12 * cos(1.0), 0.12 * sin(1.0), 1, 0, 0.04,
                        3999, 8, 8, 0.1, px, py);
  check("like charges: line escapes to the boundary",
        fabs(px[n - 1]) > 7.9 || fabs(py[n - 1]) > 7.9, 1, 0);
}

static void test_optics(void) {
  puts("optics: Snell, TIR, Brewster, Fresnel, prism, thin lens");
  double tt, d2r = M_PI / 180.0;
  check("Snell: air -> glass 30 deg",
        pil_snell(1.0, 1.5, 30 * d2r, &tt) ? tt : -1, asin(0.5 / 1.5), 1e-12);
  check("critical angle glass -> air = 41.81 deg",
        pil_critical_angle(1.5, 1.0) / d2r, 41.8103, 1e-3);
  check("no critical angle going into denser medium",
        isnan(pil_critical_angle(1.0, 1.5)), 1, 0);
  check("TIR beyond the critical angle", pil_snell(1.5, 1.0, 45 * d2r, &tt), 0,
        0);
  double bw = pil_brewster_angle(1.0, 1.5);
  check("Brewster: R_p = 0", pil_fresnel_rp(1.0, 1.5, bw), 0.0, 1e-12);
  check("Brewster: reflected and refracted rays are perpendicular",
        bw + asin(sin(bw) / 1.5), M_PI / 2, 1e-12);
  check("normal incidence R = ((n1 - n2)/(n1 + n2))^2",
        pil_fresnel_rs(1.0, 1.5, 0.0), 0.04, 1e-12);
  check("R_s = R_p at normal incidence", pil_fresnel_rp(1.0, 1.5, 0.0),
        pil_fresnel_rs(1.0, 1.5, 0.0), 1e-12);
  check("R = 1 in total internal reflection",
        pil_fresnel_rs(1.5, 1.0, 50 * d2r), 1.0, 0);
  // Vector form agrees with the angle form
  double ti = 25 * d2r;
  double d[2] = {sin(ti), -cos(ti)};
  double nrm[2] = {0, 1};
  double out[2];
  check("vector refraction gives Snell's angle",
        pil_refract(d, nrm, 1.0, 1.5, out) ? atan2(out[0], -out[1]) : -1,
        asin(sin(ti) / 1.5), 1e-12);
  check("refracted direction is a unit vector", hypot(out[0], out[1]), 1.0,
        1e-12);
  double d3[2] = {sin(50 * d2r), cos(50 * d2r)};
  double n3[2] = {0, -1};
  check("vector form: TIR returns reflection",
        pil_refract(d3, n3, 1.5, 1.0, out), 0, 0);
  check("reflected ray keeps the angle (x)", out[0], d3[0], 1e-12);
  check("reflected ray flips the normal component (y)", out[1], -d3[1], 1e-12);
  // Prism: scan the incidence angle for the minimum deviation
  double apex = 60 * d2r;
  double best = 1e9;
  double tbest = 0;
  for (int i = 0; i < 8000; i++) {
    double a = (10.0 + 0.01 * i) * d2r;
    double dv = pil_prism_deviation(apex, 1.5, a);
    if (!isnan(dv) && dv < best) {
      best = dv;
      tbest = a;
    }
  }

  check("prism min deviation matches 2 asin(n sin(A/2)) - A", best,
        pil_prism_min_deviation(apex, 1.5), 1e-6);
  check("minimum deviation at symmetric passage", asin(1.5 * sin(0.5 * apex)),
        tbest, 1e-3);
  check("thin lens: object at 2f images at 2f", pil_lens_image(0.1, 0.2), 0.2,
        1e-12);
  check("thin lens: magnification -1 at 2f", pil_lens_magnification(0.1, 0.2),
        -1.0, 1e-12);
  check("object inside f: virtual, upright, magnified",
        pil_lens_image(0.1, 0.05) < 0 && pil_lens_magnification(0.1, 0.05) > 1,
        1, 0);
}

static void test_lens_rays(void) {
  puts("thin lens: principal rays meet at the image");
  double f = 1.0;
  double h = 0.55;
  double dos[5] = {
      3.4, 2.0, 1.5, 0.8,
      0.4}; // Real inverted, 2f, between f and 2f, virtual, virtual
  for (int i = 0; i < 5; i++) {
    double d = dos[i];
    double s[3];
    double y[3];
    double s2[3];
    pil_lens_principal_slopes(f, d, h, s);
    for (int k = 0; k < 3; k++) {
      y[k] = h + s[k] * d; // Height at lens
      s2[k] = pil_lens_ray_slope(f, y[k], s[k]);
    }

    // Intersection of rays 0 and 2 after lens (lines y = y_k + s2_k x)
    double xi = (y[2] - y[0]) / (s2[0] - s2[2]);
    double yi = y[0] + s2[0] * xi;
    check("rays 0 and 2 meet at x = d_i", xi, pil_lens_image(f, d), 1e-9);
    check("ray 1 (through the centre) passes through the same point",
          y[1] + s2[1] * xi, yi, 1e-9);
    check("image height = m h", yi, pil_lens_magnification(f, d) * h, 1e-9);
  }

  double s[3];
  pil_lens_principal_slopes(f, 2.0, 0.5, s);
  check("ray 0 is parallel to the axis", s[0], 0.0, 0);
  check("parallel ray goes through the far focus",
        pil_lens_ray_slope(f, 0.5, 0.0), -0.5, 1e-12);
  check("central ray is undeflected", pil_lens_ray_slope(f, 0.0, 0.3), 0.3,
        1e-12);
}

static void test_prism_glass(void) {
  puts("prism: glass ray trace matches the closed-form deviation");
  static const double tri[3][2] = {{-1.5, 0.0}, {1.5, 0.0}, {0.0, 2.598076}};
  double apex = M_PI / 3.0;
  double n_v = pil_cauchy_index(1.60, 0.010, 0.4),
         n_r = pil_cauchy_index(1.60, 0.010, 0.7);
  check("Cauchy n(400 nm)", n_v, 1.6625, 1e-9);
  check("violet bends more than red", n_v > n_r ? 1.0 : 0.0, 1.0, 0.0);
  for (int k = 0; k < 6; k++) {
    double ti = (40.0 + 7.0 * k) * M_PI / 180.0;
    double n = k % 2 ? n_v : n_r;
    double ang = -M_PI / 6.0 + ti;
    double d[2] = {cos(ang), sin(ang)}, out[2];
    double pts[8][2];
    double pe[2] = {-0.75, 1.299038};
    double st[2] = {pe[0] - 4 * d[0], pe[1] - 4 * d[1]};
    int np = pil_glass_trace(tri, 3, n, st, d, 5.0, pts, 8, out);
    check("trace has 4 points", np, 4, 0.0);
    if (np == 4) {
      double dev = acos(d[0] * out[0] + d[1] * out[1]);
      check("deviation = closed form", dev, pil_prism_deviation(apex, n, ti),
            1e-5);
    }
  }
  double r0[3];
  double r1[3];
  pil_wavelength_rgb(400.0, r0);
  pil_wavelength_rgb(700.0, r1);
  check("violet is blue-dominant", r0[2] > r0[0] ? 1.0 : 0.0, 1.0, 0.0);
  check("red is red only", r1[0] > r1[1] && r1[2] == 0.0 ? 1.0 : 0.0, 1.0, 0.0);
}

static void test_timer555(void) {
  puts("555 astable: simulated frequency and duty match the formulas");
  double r1 = 1000.0;
  double r2 = 2200.0;
  double c = 1e-6;
  pil_555 t;
  pil_555_init(&t, r1, r2, c, 9.0);
  double h = 1e-5;
  double tm = 0.0;
  double rise = -1.0;
  double first = -1.0;
  double hi = 0.0;
  double vmin = 9.0;
  double vmax = 0.0;
  int n = 0;
  int prev = t.out;
  while (tm < 0.1) {
    pil_555_step(&t, h);
    tm += h;
    if (t.out) {
      hi += h;
    }
    if (tm > 0.01) {
      if (t.vc < vmin) {
        vmin = t.vc;
      }
      if (t.vc > vmax) {
        vmax = t.vc;
      }
    }
    if (t.out && !prev) {
      if (first < 0) {
        first = tm;
      }
      rise = tm;
      n++;
    }
    prev = t.out;
  }

  double f = (n - 1) / (rise - first);
  check("frequency", f / pil_555_frequency(r1, r2, c), 1.0, 5e-3);
  check("duty", hi / tm, pil_555_duty(r1, r2), 1e-2);
  check("v_C stays between 1/3 and 2/3 Vcc (min)", vmin, 3.0, 1e-2);
  check("v_C stays between 1/3 and 2/3 Vcc (max)", vmax, 6.0, 1e-2);
  check("duty formula = ton / T",
        pil_555_high_time(r1, r2, c) * pil_555_frequency(r1, r2, c),
        pil_555_duty(r1, r2), 1e-12);
}

static void test_equipotential(void) {
  puts("equipotentials: stay on V = const, close around one charge, V = 0 "
       "plane of a dipole");
  pil_charges c;
  c.n = 1;
  c.q[0] = 2.0;
  c.x[0][0] = 0.3;
  c.x[0][1] = -0.2;
  static double px[8004], py[8004];
  int n = pil_es_equipotential(&c, 1.3, -0.2, 0.05, 4000, 10, 10, 0.05, px, py);
  double rmax = 0;
  double rmin = 1e9;
  for (int i = 0; i < n; i++) {
    double r = hypot(px[i] - 0.3, py[i] + 0.2);
    if (r > rmax) {
      rmax = r;
    }
    if (r < rmin) {
      rmin = r;
    }
  }

  check("circle radius q/V", 0.5 * (rmax + rmin), 1.0, 1e-4);
  check("circle is round", rmax - rmin, 0.0, 1e-4);
  check("loop closes on its start", hypot(px[n - 1] - px[0], py[n - 1] - py[0]),
        0.0, 1e-9);
  c.n = 2;
  c.q[0] = 1.0;
  c.q[1] = -1.0;
  c.x[0][0] = -1.2;
  c.x[0][1] = 0.0;
  c.x[1][0] = 1.2;
  c.x[1][1] = 0.0;
  n = pil_es_equipotential(&c, 0.0, 0.0, 0.05, 400, 5, 3, 0.05, px, py);
  double ylo = 1e9;
  double yhi = -1e9;
  for (int i = 0; i < n; i++) {
    if (py[i] < ylo) {
      ylo = py[i];
    }
    if (py[i] > yhi) {
      yhi = py[i];
    }
  }
  check("x = 0 line spans the box both ways",
        ylo < -2.9 && yhi > 2.9 ? 1.0 : 0.0, 1.0, 0.0);
  double worst = 0;
  double vv = 0;
  for (int i = 0; i < n; i++) {
    if (fabs(px[i]) > worst) {
      worst = fabs(px[i]);
    }
    double v = fabs(pil_es_potential(&c, px[i], py[i]));
    if (v > vv) {
      vv = v;
    }
  }

  check("V = 0 curve is the line x = 0", worst, 0.0, 1e-9);
  check("potential stays 0 along it", vv, 0.0, 1e-9);
}

static void test_solenoid(void) {
  puts(
      "solenoid: summed loops vs the exact sheet formula, long-solenoid limit");
  double bx;
  double bz;
  for (int k = 0; k < 3; k++) {
    double z = 0.0 + 1.1 * k;
    pil_solenoid(1.0, 6.0, 24, 0.0, z, 200, &bx, &bz);
    check("on-axis sum vs sheet formula",
          bz / pil_solenoid_axis(1.0, 6.0, 4.0, z), 1.0, 3e-3);
    check("axial field only on the axis", bx, 0.0, 1e-9);
  }

  pil_solenoid(1.0, 20.0, 80, 0.0, 0.0, 200, &bx, &bz);
  double bin = bz;
  check("long solenoid centre -> 4 pi n", bin / pil_solenoid_ideal(4.0), 1.0,
        6e-3);
  pil_solenoid(1.0, 20.0, 80, 3.0, 0.0, 200, &bx, &bz);
  check("outside field is < 1 % of inside", fabs(bz) / bin < 0.01 ? 1.0 : 0.0,
        1.0, 0.0);
  check("sheet formula centre = n 4 pi L / sqrt(L^2 + 4 r^2)",
        pil_solenoid_axis(1.0, 6.0, 4.0, 0.0) / pil_solenoid_ideal(4.0),
        6.0 / sqrt(36.0 + 4.0), 1e-12);
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
  test_top();
  test_rigid3();
  test_rolling_on_plank();
  test_trace();
  test_central_force();
  test_lorentz_drag();
  test_collision();
  test_modes();
  test_bead_chain();
  test_roll_race();
  test_bridge_rectifier();
  test_resonance();
  test_devices();
  test_hoop();
  test_double_pendulum();
  test_kirchhoff();
  test_doppler();
  test_spherical();
  test_mosfet();
  test_rlc_response();
  test_boris();
  test_dipole();
  test_foucault();
  test_chaos();
  test_helmholtz();
  test_hohmann();
  test_nbody();
  test_emwave();
  test_rutherford();
  test_electrostatic();
  test_optics();
  test_lens_rays();
  test_prism_glass();
  test_timer555();
  test_equipotential();
  test_solenoid();
  printf("\n%s (%d failure%s)\n", fails ? "FAILED" : "ALL PASSED", fails,
         fails == 1 ? "" : "s");
  return fails ? 1 : 0;
}
