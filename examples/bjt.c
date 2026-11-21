/*
 * bjt.c - Bipolar Junction Transistor; common-emitter amplifier
 * R_B sets bias (I_B ~ 20 uA, I_C = beta I_B ~ 2 mA, V_C ~ 7.6 V);
 * coupling capacitor passes a 5 mV, 1 kHz signal onto base and collector swings
 * opposite way with gain ~ -g_m R_C
 * Exponential transistor curve makes output slightly asymmetric
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.001
#define PIL_EXAMPLE_IMAGE_TIME 0.0023
#include "example_main.h"
#include <math.h>
#include <pil/pil_circuit.h>

#define VCC 12.0

// Collector swing about its bias point, in volts
static double vc_ac(const pil_world *w, void *user) {
  (void)w;

  double v = pil_cdevice_voltage((const pil_cdevice *)user, 0, 2);

  return v == 0.0 ? 0.0 : v - 7.6; // 0 before first solve
}

static double vin_x100(const pil_world *w, void *user) {
  (void)w;

  return 100.0 * pil_celement_voltage((const pil_celement *)user);
}

static double ic_ma(const pil_world *w, void *user) {
  (void)w;

  return 1e3 * pil_cdevice_current((const pil_cdevice *)user, 0);
}

int main(int argc, char **argv) {
  pil_world *w = pil_world_new();
  pil_circuit *c = pil_add_circuit(w);

  pil_node a = pil_circuit_node(c, 0.0, 1.0), g = pil_circuit_node(c, 0.0, 0.0),
           nb = pil_circuit_node(c, 1.6, 1.0),
           rt = pil_circuit_node(c, 1.6, 3.2),
           vcc = pil_circuit_node(c, 3.2, 3.2),
           col = pil_circuit_node(c, 3.2, 2.0),
           bas = pil_circuit_node(c, 2.4, 1.0),
           em = pil_circuit_node(c, 3.2, 0.0),
           vt = pil_circuit_node(c, 4.8, 3.2),
           gb = pil_circuit_node(c, 4.8, 0.0);
  pil_circuit_ground(c, g);

  pil_celement *src = pil_vsource_ac(c, a, g, 5e-3, 1000.0, 0.0, "v_{in}");
  pil_capacitor(c, a, nb, 10e-6, -0.677, "C_c");
  pil_resistor(c, rt, nb, 560e3, "R_B");
  pil_wire(c, rt, vcc);
  pil_wire(c, nb, bas);
  pil_resistor(c, vcc, col, 2.2e3, "R_C");
  pil_cdevice *q = pil_npn(c, col, bas, em, 100.0, "Q");
  pil_wire(c, g, em);
  pil_wire(c, em, gb);
  pil_vsource_dc(c, vt, gb, VCC, "V_{CC}");
  pil_wire(c, vcc, vt);
  pil_celement_set_label_at(src, -34, 0);
  pil_circuit_set_step(c, 2e-6);

  pil_plot *p = pil_add_plot(w, 0.0, 0.56, 0.64, 1.0, "t [s]", "[V]", 0.004);
  pil_plot_add_signal(p, "100 v_{in}", vin_x100, src);
  pil_plot_add_signal(p, "v_C - V_Q", vc_ac, q);
  pil_plot *r =
      pil_add_plot(w, 0.68, 0.56, 0.97, 1.0, "t [s]", "I_C [mA]", 0.004);
  pil_plot_add_signal(r, "I_C", ic_ma, q);
  pil_plots_record(w, 1e-5, 2e-4);

  pil_view view;
  pil_view_fit_region(&view, -0.8, -0.45, 5.6, 3.6, 1100, 800, 0.0, 0.0, 1.0,
                      0.55, 18);
  int rc =
      example_main(w, &view, "Common-emitter amplifier", argc, argv, "bjt.png");

  pil_world_free(w);

  return rc;
}
