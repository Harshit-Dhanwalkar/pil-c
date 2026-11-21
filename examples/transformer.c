/*
 * transformer.c - step-up transformer, n = 2, with a resistive load
 * Two inductors coupled with pil_mutual (k = 0.998); dots mark ends whose
 * voltages are in phase
 * Secondary voltage is n times primary and secondary current is 1 / n of
 * primary's (apart from small magnetising current)
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.0033
#define PIL_EXAMPLE_IMAGE_TIME 0.0052
#include "example_main.h"
#include <math.h>
#include <pil/pil_circuit.h>

static double volts_of(const pil_world *w, void *user) {
  (void)w;

  return pil_celement_voltage((const pil_celement *)user);
}

int main(int argc, char **argv) {
  pil_world *w = pil_world_new();
  pil_circuit *c = pil_add_circuit(w);

  pil_node sa = pil_circuit_node(c, 0.0, 1.8),
           sb = pil_circuit_node(c, 0.0, 0.0),
           p1 = pil_circuit_node(c, 1.6, 1.8),
           p2 = pil_circuit_node(c, 1.6, 0.0),
           s1 = pil_circuit_node(c, 3.0, 1.8),
           s2 = pil_circuit_node(c, 3.0, 0.0),
           l1 = pil_circuit_node(c, 4.4, 1.8),
           l2 = pil_circuit_node(c, 4.4, 0.0);
  pil_circuit_ground(c, sb);

  pil_celement *src = pil_vsource_ac(c, sa, sb, 10.0, 400.0, 0.0, "v_1");
  pil_celement *prim = pil_inductor(c, p1, p2, 1.0, 0.0, "L_1");
  pil_celement *sec = pil_inductor(c, s1, s2, 4.0, 0.0, "L_2");
  pil_mutual(prim, sec, 0.998);
  pil_wire(c, sa, p1);
  pil_wire(c, sb, p2);
  pil_wire(c, s1, l1);
  pil_wire(c, s2, l2);
  pil_celement *load = pil_resistor(c, l1, l2, 2000.0, "R");
  pil_celement_set_label_at(src, -34, 0);
  pil_celement_set_label_at(prim, -34, 0);
  pil_celement_set_label_at(sec, 34, 0);
  pil_circuit_set_step(c, 5e-6);

  pil_plot *p = pil_add_plot(w, 0.0, 0.58, 1.0, 1.0, "t [s]", "v [V]", 0.02);
  pil_plot_add_signal(p, "v_1", volts_of, src);
  pil_plot_add_signal(p, "v_2", volts_of, load);
  pil_plots_record(w, 1e-4, 5e-4);

  pil_view view;
  pil_view_fit_region(&view, -0.8, -0.45, 5.0, 2.1, 900, 760, 0.0, 0.0, 1.0,
                      0.57, 18);
  int rc = example_main(w, &view, "Transformer", argc, argv, "transformer.png");

  pil_world_free(w);

  return rc;
}
