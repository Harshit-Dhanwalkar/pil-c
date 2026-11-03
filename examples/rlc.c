/*
 * rlc.c - a driven series RLC circuit, with voltages plotted beneath it
 * Source v(t) -> R -> L -> C, all in one loop; panel shows source, resistor and
 * capacitor voltage
 *Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.0033 // 1 s of video = 3.3 ms of circuit
#include "example_main.h"
#include <math.h>
#include <pil/pil_circuit.h>

static double volts_of(const pil_world *w, void *user) {
  (void)w;

  return pil_celement_voltage((const pil_celement *)user);
}

int main(int argc, char **argv) {
  // f0 = 503 Hz, Q = 1.6
  const double R = 20.0;
  const double L = 10e-3;
  const double C = 10e-6;
  const double f = 400.0;
  const double V = 10.0;

  pil_world *w = pil_world_new();
  pil_circuit *c = pil_add_circuit(w);

  pil_node tl = pil_circuit_node(c, 0.0, 1.6),
           tm = pil_circuit_node(c, 2.0, 1.6),
           tr = pil_circuit_node(c, 4.0, 1.6),
           br = pil_circuit_node(c, 4.0, 0.0),
           bl = pil_circuit_node(c, 0.0, 0.0);
  pil_circuit_ground(c, bl);

  pil_celement *src = pil_vsource_ac(c, tl, bl, V, f, 0.0, "v(t)");
  pil_celement *res = pil_resistor(c, tl, tm, R, "R");
  pil_inductor(c, tm, tr, L, 0.0, "L");
  pil_celement *cap = pil_capacitor(c, tr, br, C, 0.0, "C");
  pil_wire(c, bl, br);
  pil_celement_set_label_at(src, -34, 0);

  pil_add_arrow(w, pil_at(1.0, 1.6), 46, 0, "i");

  pil_plot *p = pil_add_plot(w, 0.0, 0.58, 1.0, 1.0, "t [s]", "v [V]", 0.02);
  pil_plot_add_signal(p, "v", volts_of, src);
  pil_plot_add_signal(p, "v_R", volts_of, res);
  pil_plot_add_signal(p, "v_C", volts_of, cap);
  pil_plots_record(w, 1e-4, 1e-3);

  pil_view view;
  pil_view_fit_region(&view, -0.7, -0.45, 4.9, 2.0, 900, 760, 0.0, 0.0, 1.0,
                      0.57, 18);
  int rc = example_main(w, &view, "Series RLC circuit", argc, argv, "rlc.png");

  pil_world_free(w);

  return rc;
}
