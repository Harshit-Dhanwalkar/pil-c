/*
 * zener.c - zener shunt regulator
 * A 12 V supply with 3 V of 100 Hz ripple feeds a 5.1 V zener through
 * R_s = 470 ohm; in breakdown zener clamps load voltage so ripple is reduced by
 * r_z / (R_s + r_z), here roughly a factor 50
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.004
#define PIL_EXAMPLE_IMAGE_TIME 0.0065
#include "example_main.h"
#include <math.h>
#include <pil/pil_circuit.h>

static double supply(double t, void *u) {
  (void)u;

  return 12.0 + 3.0 * sin(2.0 * M_PI * 100.0 * t);
}

static double volts_of(const pil_world *w, void *user) {
  (void)w;

  return pil_celement_voltage((const pil_celement *)user);
}

static double amps_of(const pil_world *w, void *user) {
  (void)w;

  return 1e3 * pil_celement_current((const pil_celement *)user);
}

// Zener is wired anode -> cathode upwards, so its reverse current is negative
static double amps_rev(const pil_world *w, void *user) {
  return -amps_of(w, user);
}

int main(int argc, char **argv) {
  pil_world *w = pil_world_new();
  pil_circuit *c = pil_add_circuit(w);

  pil_node a = pil_circuit_node(c, 0.0, 1.8), g = pil_circuit_node(c, 0.0, 0.0),
           z = pil_circuit_node(c, 2.2, 1.8),
           zg = pil_circuit_node(c, 2.2, 0.0),
           l = pil_circuit_node(c, 3.8, 1.8),
           lg = pil_circuit_node(c, 3.8, 0.0);
  pil_circuit_ground(c, g);

  pil_celement *src = pil_vsource_fn(c, a, g, supply, NULL, "v_s");
  pil_celement *rs = pil_resistor(c, a, z, 470.0, "R_s");
  pil_celement *zd = pil_zener(c, zg, z, 5.1, "5.1 V");
  pil_wire(c, z, l);
  pil_celement *load = pil_resistor(c, l, lg, 2000.0, "R_L");
  pil_wire(c, g, zg);
  pil_wire(c, zg, lg);
  pil_celement_set_label_at(src, -34, 0);
  pil_celement_set_label_at(zd, 40, 0);
  pil_circuit_set_step(c, 5e-6);

  pil_plot *p = pil_add_plot(w, 0.0, 0.58, 0.55, 1.0, "t [s]", "v [V]", 0.03);
  pil_plot_add_signal(p, "v_s", volts_of, src);
  pil_plot_add_signal(p, "v_L", volts_of, load);
  pil_plot *q = pil_add_plot(w, 0.62, 0.58, 0.97, 1.0, "t [s]", "i [mA]", 0.03);
  pil_plot_add_signal(q, "i_R", amps_of, rs);
  pil_plot_add_signal(q, "i_Z", amps_rev, zd);
  pil_plots_record(w, 1e-4, 5e-4);

  pil_view view;
  pil_view_fit_region(&view, -0.8, -0.45, 4.5, 2.1, 1000, 760, 0.0, 0.0, 1.0,
                      0.57, 18);
  int rc = example_main(w, &view, "Zener regulator", argc, argv, "zener.png");

  pil_world_free(w);

  return rc;
}
