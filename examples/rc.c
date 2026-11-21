/*
 * rc.c - RC charging and discharging driven by a square wave
 * v_C = V (1 - exp(-t / RC)) while source is on, V exp(-(t - T/2) / RC) after
 * it is switched to 0 V
 * Time constant tau = RC = 1 ms; plot marks v_C = 63 % of V at one tau with
 * guide curves of ideal law
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.0008
#define PIL_EXAMPLE_IMAGE_TIME 0.0055
#include "example_main.h"
#include <math.h>
#include <pil/pil_circuit.h>

#define V0 10.0
#define R_ 1000.0
#define C_ 1e-6
#define HALF 4e-3 // Source on for 4 ms, off for 4 ms

static double square(double t, void *u) {
  (void)u;

  return fmod(t, 2.0 * HALF) < HALF ? V0 : 0.0;
}

static double volts_of(const pil_world *w, void *user) {
  (void)w;

  return pil_celement_voltage((const pil_celement *)user);
}

static double amps_of(const pil_world *w, void *user) {
  (void)w;

  return 1e3 * pil_celement_current((const pil_celement *)user); // mA
}

int main(int argc, char **argv) {
  pil_world *w = pil_world_new();
  pil_circuit *c = pil_add_circuit(w);
  pil_node tl = pil_circuit_node(c, 0.0, 1.6),
           tr = pil_circuit_node(c, 3.0, 1.6),
           br = pil_circuit_node(c, 3.0, 0.0),
           bl = pil_circuit_node(c, 0.0, 0.0);
  pil_circuit_ground(c, bl);

  pil_celement *src = pil_vsource_fn(c, tl, bl, square, NULL, "v(t)");
  pil_celement *res = pil_resistor(c, tl, tr, R_, "R");
  pil_celement *cap = pil_capacitor(c, tr, br, C_, 0.0, "C");
  pil_wire(c, bl, br);
  pil_celement_set_label_at(src, -34, 0);
  pil_add_arrow(w, pil_at(1.5, 1.6), 46, 0, "i");

  pil_plot *p = pil_add_plot(w, 0.0, 0.58, 0.6, 1.0, "t [s]", "v [V]", 0.008);
  pil_plot_add_signal(p, "v", volts_of, src);
  pil_plot_add_signal(p, "v_C", volts_of, cap);
  pil_plot *q = pil_add_plot(w, 0.62, 0.58, 1.0, 1.0, "t [s]", "i [mA]", 0.008);
  pil_plot_add_signal(q, "i", amps_of, res);
  pil_plots_record(w, 1e-5, 1e-4);

  pil_view view;
  pil_view_fit_region(&view, -0.7, -0.45, 3.7, 2.0, 1000, 760, 0.0, 0.0, 1.0,
                      0.57, 18);
  int rc = example_main(w, &view, "RC charging", argc, argv, "rc.png");

  pil_world_free(w);

  return rc;
}
