/*
 * cmos.c - CMOS inverter and its transfer characteristic
 * An nmos pulls output down, a pmos pulls it up; exactly one conducts once
 * input is far from V_dd / 2, so output sits at a rail and gate draws almost no
 * power
 * A triangle wave on input sweeps gate through switching threshold: output
 * against input (right) is sharp transfer curve and capacitor C_L is charged /
 * discharged by transistors between rails
 */

#define PIL_EXAMPLE_TIME_SCALE 0.001
#define PIL_EXAMPLE_IMAGE_TIME 0.0023
#include "example_main.h"
#include <math.h>
#include <pil/pil_circuit.h>

#define VDD 5.0
#define PERIOD 4e-3

static double tri(double t, void *u) {
  (void)u;
  double ph = fmod(t, PERIOD) / PERIOD;

  return VDD * (ph < 0.5 ? 2.0 * ph : 2.0 - 2.0 * ph);
}

static double node_v(const pil_world *w, void *u) {
  (void)w;
  const pil_celement *e = u;

  return pil_celement_voltage(e);
}

static double out_v(const pil_world *w, void *u) {
  pil_circuit *c = ((void **)u)[0];
  pil_node n = *(pil_node *)((void **)u)[1];

  return pil_world_time(w) == 0.0
             ? VDD
             : pil_circuit_node_voltage(c, n); // C_L starts charged
}

static double idd_ma(const pil_world *w, void *u) {
  (void)w;

  return -1e3 *
         pil_cdevice_current((const pil_cdevice *)u, 2); // pmos source current
}

int main(int argc, char **argv) {
  pil_world *w = pil_world_new();
  pil_circuit *c = pil_add_circuit(w);

  pil_node sa = pil_circuit_node(c, 0.0, 1.8),
           g = pil_circuit_node(c, 0.0, 0.0),
           in = pil_circuit_node(c, 1.0, 1.8),
           gpa = pil_circuit_node(c, 1.0, 2.7),
           gna = pil_circuit_node(c, 1.0, 0.9),
           gp = pil_circuit_node(c, 1.9, 2.7),
           gn = pil_circuit_node(c, 1.9, 0.9),
           vdd = pil_circuit_node(c, 2.8, 3.6),
           out = pil_circuit_node(c, 2.8, 1.8),
           outr = pil_circuit_node(c, 4.4, 1.8),
           vt = pil_circuit_node(c, 5.8, 3.6),
           gr = pil_circuit_node(c, 5.8, 0.0),
           sn = pil_circuit_node(c, 2.8, 0.0),
           gc = pil_circuit_node(c, 4.4, 0.0);
  pil_circuit_ground(c, g);

  pil_celement *vin = pil_vsource_fn(c, sa, g, tri, NULL, "v_{in}");
  pil_wire(c, sa, in);
  pil_wire(c, in, gpa);
  pil_wire(c, in, gna);
  pil_wire(c, gpa, gp);
  pil_wire(c, gna, gn);
  pil_cdevice *mp = pil_pmos(c, out, gp, vdd, 2e-3, 1.0, "M_p");
  pil_cdevice *mn = pil_nmos(c, out, gn, sn, 2e-3, 1.0, "M_n");
  (void)mn;
  pil_wire(c, g, sn);
  pil_wire(c, out, outr);
  pil_capacitor(c, outr, gc, 10e-9, VDD, "C_L");
  pil_wire(c, sn, gc);
  pil_wire(c, gc, gr);
  pil_vsource_dc(c, vt, gr, VDD, "V_{dd}");
  pil_wire(c, vdd, vt);
  pil_celement_set_label_at(vin, -34, 0);
  pil_circuit_set_step(c, 2e-6);

  static pil_node outn;
  outn = out;
  static void *ref[2];
  ref[0] = c;
  ref[1] = &outn;

  pil_plot *p = pil_add_plot(w, 0.0, 0.60, 0.42, 1.0, "t [s]", "v [V]", 0.008);
  pil_plot_add_signal(p, "v_{in}", node_v, vin);
  pil_plot_add_signal(p, "v_{out}", out_v, ref);
  pil_plot *q = pil_add_plot(w, 0.49, 0.60, 0.70, 1.0, "v_{in} [V]",
                             "v_{out} [V]", 0.008);
  pil_plot_add_signal(q, "", node_v, vin);
  pil_plot_add_signal(q, "", out_v, ref);
  pil_plot_set_xy(q, 0);
  pil_plot *r =
      pil_add_plot(w, 0.77, 0.60, 0.98, 1.0, "t [s]", "I_{dd} [mA]", 0.008);
  pil_plot_add_signal(r, "", idd_ma, mp);
  pil_plots_record(w, 2e-5, 2e-4);

  pil_view view;
  pil_view_fit_region(&view, -0.7, -0.45, 6.4, 4.0, 1300, 760, 0.0, 0.0, 1.0,
                      0.58, 18);
  int rc = example_main(w, &view, "CMOS inverter", argc, argv, "cmos.png");

  pil_world_free(w);

  return rc;
}
