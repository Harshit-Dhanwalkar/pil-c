/*
 * lc.c - LC tank: energy sloshes between capacitor and inductor
 * C starts charged to 10 V; with a little resistance oscillation at 
 * f0 = 1 / (2 pi sqrt(LC)) = 503 Hz decays
 * Energy panel: E_C and E_L trade places, their sum falls only through R
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.0015
#define PIL_EXAMPLE_IMAGE_TIME 0.0021
#include "example_main.h"
#include <math.h>
#include <pil/pil_circuit.h>

#define L_ 10e-3
#define C_ 10e-6

static double volts_of(const pil_world *w, void *user) {
  (void)w;

  return pil_celement_voltage((const pil_celement *)user);
}

static double amps_of(const pil_world *w, void *user) {
  (void)w;

  return 10.0 * pil_celement_current((const pil_celement *)user); // 0.1 A units
}

static double eC(const pil_world *w, void *user) {
  (void)w;
  double v = pil_celement_voltage((const pil_celement *)user);

  return 1e3 * 0.5 * C_ * v * v; // mJ
}

static double eL(const pil_world *w, void *user) {
  (void)w;
  double i = pil_celement_current((const pil_celement *)user);

  return 1e3 * 0.5 * L_ * i * i;
}

int main(int argc, char **argv) {
  pil_world *w = pil_world_new();
  pil_circuit *c = pil_add_circuit(w);
  pil_node tl = pil_circuit_node(c, 0.0, 1.6),
           tm = pil_circuit_node(c, 1.6, 1.6),
           tr = pil_circuit_node(c, 3.2, 1.6),
           br = pil_circuit_node(c, 3.2, 0.0),
           bl = pil_circuit_node(c, 0.0, 0.0);
  pil_circuit_ground(c, bl);

  pil_celement *cap = pil_capacitor(c, tl, bl, C_, 10.0, "C");
  pil_resistor(c, tl, tm, 2.0, "R");
  pil_celement *ind = pil_inductor(c, tm, tr, L_, 0.0, "L");
  pil_wire(c, tr, br);
  pil_wire(c, br, bl);

  pil_add_arrow(w, pil_at(2.4, 1.6), 46, 0, "i");

  pil_plot *p = pil_add_plot(w, 0.0, 0.58, 0.6, 1.0, "t [s]", "", 0.004);
  pil_plot_add_signal(p, "v_C [V]", volts_of, cap);
  pil_plot_add_signal(p, "10 i [A]", amps_of, ind);
  pil_plot *q = pil_add_plot(w, 0.62, 0.58, 1.0, 1.0, "t [s]", "E [mJ]", 0.004);
  pil_plot_add_signal(q, "E_C", eC, cap);
  pil_plot_add_signal(q, "E_L", eL, ind);
  pil_plots_record(w, 1e-5, 1e-4);

  pil_view view;
  pil_view_fit_region(&view, -0.7, -0.45, 3.9, 2.0, 1000, 760, 0.0, 0.0, 1.0,
                      0.57, 18);
  int rc = example_main(w, &view, "LC oscillation", argc, argv, "lc.png");

  pil_world_free(w);

  return rc;
}
