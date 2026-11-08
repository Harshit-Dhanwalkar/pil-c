/*
 * motor.c - a DC motor spins up a flywheel: circuit and mechanism in one world
 * Battery -> switch -> motor; shaft turns a flywheel J against viscous friction
 * Motor torque is k_t i and it feeds back k_t omega as back-EMF, which plots
 * show as inrush current decaying while speed climbs
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.06 // 1 s of video = 60 ms of motion
#define PIL_EXAMPLE_IMAGE_TIME 0.03
#include "example_main.h"
#include <math.h>
#include <pil/pil_circuit.h>

static double omega_of(const pil_world *w, void *user) {
  (void)w;

  return ((const pil_body *)user)->omega;
}

static double amps_of(const pil_world *w, void *user) {
  (void)w;

  return pil_celement_current((const pil_celement *)user);
}

int main(int argc, char **argv) {
  // Motor and battery/
  const double kt = 0.05;
  const double R = 1.0;
  const double L = 1e-3;
  const double V = 5.0;
  // Flywheel, friction
  const double J = 2e-4;
  const double c = 1e-3;

  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);

  // Flywheel: a disc on a fixed axle, with viscous drag
  pil_vec2 axle = {4.0, 0.6};
  pil_body *wheel = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(0.6),
                                                     .mass = 1.0,
                                                     .inertia = J,
                                                     .pos = axle,
                                                     .flags = PIL_SPOKE,
                                                     .label = "J",
                                                     .label_dx = 40,
                                                     .label_dy = 40});
  pil_add_revolute(w, pil_on(wheel, 0, 0), pil_at(axle.x, axle.y));
  pil_add_torsion_spring(w, wheel, NULL, 0.0, 0.0, c);

  // Circuit
  pil_circuit *ck = pil_add_circuit(w);
  pil_node tl = pil_circuit_node(ck, 0.0, 1.2),
           tr = pil_circuit_node(ck, 2.0, 1.2),
           br = pil_circuit_node(ck, 2.0, 0.0),
           bl = pil_circuit_node(ck, 0.0, 0.0);
  pil_circuit_ground(ck, bl);
  pil_vsource_dc(ck, tl, bl, V, "V");
  pil_switch(ck, tl, tr, 0.05, 0, "S");
  pil_celement *motor = pil_motor(ck, tr, br, wheel, NULL, kt, R, L, NULL);
  pil_wire(ck, bl, br);

  // shaft from motor to flywheel, and its direction of rotation
  pil_add_line(w, pil_at(2.0 + 0.22, 0.6), pil_at(axle.x, axle.y), 0);
  pil_add_arrow(w, pil_at(axle.x - 0.35, axle.y + 0.8), 30, 18, "\xcf\x89");

  pil_plot *p1 =
      pil_add_plot(w, 0.0, 0.52, 1.0, 0.76, "", "\xcf\x89 [rad/s]", 0.4);
  pil_plot_add_signal(p1, "", omega_of, wheel);
  pil_plot *p2 = pil_add_plot(w, 0.0, 0.76, 1.0, 1.0, "t [s]", "i [A]", 0.4);
  pil_plot_add_signal(p2, "", amps_of, motor);
  pil_plots_record(w, 2e-3, 1e-4);

  pil_view view;
  pil_view_fit_region(&view, -0.7, -0.4, 4.8, 1.9, 900, 800, 0.0, 0.0, 1.0, 0.5,
                      18);
  int rc =
      example_main(w, &view, "DC motor and flywheel", argc, argv, "motor.png");

  pil_world_free(w);

  return rc;
}
