/*
 * rectifier.c - half-wave rectifier with a smoothing capacitor
 * v(t) -> diode -> (load R || C); panel shows input and smoothed output
 *Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.0167 /* 1 s of video = 16.7 ms of circuit */
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

  pil_node tl = pil_circuit_node(c, 0.0, 1.5),
           ta = pil_circuit_node(c, 1.6, 1.5),
           tc = pil_circuit_node(c, 2.4, 1.5),
           tr = pil_circuit_node(c, 3.6, 1.5),
           bl = pil_circuit_node(c, 0.0, 0.0),
           bc = pil_circuit_node(c, 2.4, 0.0),
           br = pil_circuit_node(c, 3.6, 0.0);
  pil_circuit_ground(c, bl);

  pil_celement *src = pil_vsource_ac(c, tl, bl, 5.0, 50.0, 0.0, "v(t)");
  pil_diode(c, tl, ta, 1e-12, 0.026, "D");
  pil_celement *cap = pil_capacitor(c, tc, bc, 100e-6, 0.0, "C");
  pil_resistor(c, tr, br, 500.0, "R");
  pil_wire(c, ta, tc);
  pil_wire(c, tc, tr);
  pil_wire(c, bl, bc);
  pil_wire(c, bc, br);
  pil_celement_set_label_at(src, -34, 0);

  pil_plot *p = pil_add_plot(w, 0.0, 0.55, 1.0, 1.0, "t [s]", "v [V]", 0.1);
  pil_plot_add_signal(p, "v", volts_of, src);
  pil_plot_add_signal(p, "v_C", volts_of, cap);
  pil_plots_record(w, 5e-4, 1e-3);

  pil_view view;
  pil_view_fit_region(&view, -0.7, -0.45, 4.3, 1.95, 900, 760, 0.0, 0.0, 1.0,
                      0.53, 18);
  int rc = example_main(w, &view, "Half-wave rectifier", argc, argv,
                        "rectifier.png");

  pil_world_free(w);

  return rc;
}
