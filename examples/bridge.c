/*
 * bridge.c - full-wave bridge rectifier with a smoothing capacitor
 * AC source on left feeds a diode diamond; load R and filter C sit inside it
 * Capacitor is topped up twice per mains cycle, so v_load ripples at 2 f
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.2
#define PIL_EXAMPLE_IMAGE_TIME 20.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_circuit.h>

static double volts_of(const pil_world *w, void *user) {
  (void)w;

  return pil_celement_voltage((const pil_celement *)user);
}

int main(int argc, char **argv) {
  const double IS = 1e-12;
  const double NVT = 0.04;
  pil_world *w = pil_world_new();
  pil_circuit *c = pil_add_circuit(w);

  pil_node a = pil_circuit_node(c, -1.0, 1.0),
           b = pil_circuit_node(c, -1.0, -1.0),
           L = pil_circuit_node(c, 0.6, 1.0),
           Rn = pil_circuit_node(c, 4.2, 1.0),
           p = pil_circuit_node(c, 2.4, 2.4),
           n = pil_circuit_node(c, 2.4, -0.4),
           bb = pil_circuit_node(c, 5.4, -1.0),
           bt = pil_circuit_node(c, 5.4, 1.0),
           pm = pil_circuit_node(c, 2.4, 1.8),
           pa = pil_circuit_node(c, 1.9, 1.8),
           pb = pil_circuit_node(c, 2.9, 1.8),
           na = pil_circuit_node(c, 1.9, 0.6),
           nb = pil_circuit_node(c, 2.9, 0.6),
           nm = pil_circuit_node(c, 2.4, 0.6),
           bm = pil_circuit_node(c, -1.0, -1.0);
  (void)bm;
  pil_circuit_ground(c, n);

  // Dropped frequency to 0.2 Hz (period = 5s) so the arrows move slowly
  pil_celement *src = pil_vsource_ac(c, a, b, 12.0, 0.2, 0.0, "v(t)");
  pil_wire(c, a, L);
  pil_wire(c, b, bb);
  pil_wire(c, bb, bt);
  pil_wire(c, bt, Rn);
  pil_diode(c, L, p, IS, NVT, "");
  pil_diode(c, Rn, p, IS, NVT, "");
  pil_diode(c, n, L, IS, NVT, "");
  pil_diode(c, n, Rn, IS, NVT, "");
  pil_wire(c, p, pm);
  pil_wire(c, pm, pa);
  pil_wire(c, pm, pb);
  pil_celement *load = pil_resistor(c, pa, na, 500.0, "R");
  // Increased capacitor to 10,000 uF so the RC time constant matches 0.2 Hz
  pil_capacitor(c, pb, nb, 10000e-6, 0.0, "C");
  pil_wire(c, na, nm);
  pil_wire(c, nb, nm);
  pil_wire(c, nm, n);
  pil_celement_set_label_at(src, -34, 0);

  // Span increased to 10s so the plot shows two full AC cycles
  pil_plot *pl = pil_add_plot(w, 0.0, 0.58, 1.0, 1.0, "t [s]", "v [V]", 10.0);
  pil_plot_add_signal(pl, "v", volts_of, src);
  pil_plot_add_signal(pl, "v_{load}", volts_of, load);
  pil_plots_record(w, 0.05, 1e-3);

  pil_view view;
  pil_view_fit_region(&view, -1.8, -1.4, 6.0, 2.8, 1100, 800, 0.0, 0.0, 1.0,
                      0.57, 18);
  int rc = example_main(w, &view, "Bridge rectifier", argc, argv, "bridge.png");

  pil_world_free(w);

  return rc;
}
