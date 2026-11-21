/*
 * opamp.c - inverting amplifier that clips
 * Gain -R_f / R_in = -10 with +-12 V supply rails; with 1.5 V input output
 * would reach 15 V and is clipped by saturation
 * Negative feedback holds inverting input at virtual ground (node voltage ~ 0)
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.0001
#define PIL_EXAMPLE_IMAGE_TIME 0.0017
#include "example_main.h"
#include <math.h>
#include <pil/pil_circuit.h>

static double volts_of(const pil_world *w, void *user) {
  (void)w;

  return pil_celement_voltage((const pil_celement *)user);
}

static double node_of(const pil_world *w, void *user) {
  (void)w;

  return 1e3 * pil_cdevice_voltage((const pil_cdevice *)user, 1, 0); // mV
}

int main(int argc, char **argv) {
  pil_world *w = pil_world_new();
  pil_circuit *c = pil_add_circuit(w);

  pil_node a = pil_circuit_node(c, 0.0, 1.8), g = pil_circuit_node(c, 0.0, 0.0),
           sum = pil_circuit_node(c, 1.6, 1.8),
           top = pil_circuit_node(c, 1.6, 2.8),
           topr = pil_circuit_node(c, 4.0, 2.8),
           neg = pil_circuit_node(c, 2.2, 1.8),
           pos = pil_circuit_node(c, 2.2, 0.8),
           g2 = pil_circuit_node(c, 2.2, 0.0),
           out = pil_circuit_node(c, 4.0, 1.3),
           out2 = pil_circuit_node(c, 5.0, 1.3),
           gb = pil_circuit_node(c, 5.0, 0.0);
  pil_circuit_ground(c, g);

  pil_celement *src = pil_vsource_ac(c, a, g, 1.5, 1000.0, 0.0, "v_{in}");
  pil_resistor(c, a, sum, 1e3, "R_{in}");
  pil_wire(c, sum, top);
  pil_resistor(c, top, topr, 10e3, "R_f");
  pil_wire(c, topr, out);
  pil_wire(c, sum, neg);
  pil_wire(c, pos, g2);
  pil_wire(c, g, g2);
  pil_wire(c, g2, gb);
  pil_wire(c, out, out2);
  pil_celement *load = pil_resistor(c, out2, gb, 10e3, "R_L");
  pil_cdevice *amp = pil_opamp(c, pos, neg, out, 12.0, "");
  pil_celement_set_label_at(src, -34, 0);
  pil_circuit_set_step(c, 2e-6);

  pil_plot *p = pil_add_plot(w, 0.0, 0.58, 0.64, 1.0, "t [s]", "v [V]", 0.002);
  pil_plot_add_signal(p, "v_{in}", volts_of, src);
  pil_plot_add_signal(p, "v_{out}", volts_of, load);
  pil_plot *q =
      pil_add_plot(w, 0.68, 0.58, 0.97, 1.0, "t [s]", "v_- [mV]", 0.002);
  pil_plot_add_signal(q, "v_-", node_of, amp);
  pil_plots_record(w, 1e-5, 2e-4);

  pil_view view;
  pil_view_fit_region(&view, -0.8, -0.45, 5.5, 3.2, 1100, 800, 0.0, 0.0, 1.0,
                      0.57, 18);
  int rc =
      example_main(w, &view, "Inverting amplifier", argc, argv, "opamp.png");

  pil_world_free(w);

  return rc;
}
