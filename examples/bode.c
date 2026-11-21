/*
 * bode.c - frequency response of a series RLC from a swept-frequency source
 * Source frequency rises from 100 to 1000 Hz in 120 ms
 * Output taken across R is a band-pass (gain 1 at f_0 = 1 / (2 pi sqrt(LC)) =
 * 503 Hz, width f_0 / Q) and output across C a low-pass that peaks at Q ~ 3.2
 * near f_0 and falls as 1 / f^2 beyond
 * Right-hand plot shows amplitude of both against instantaneous drive frequency
 * (amplitude from quadrature pair x, x'/omega)
 */

#define PIL_EXAMPLE_TIME_SCALE 0.02
#define PIL_EXAMPLE_IMAGE_TIME 0.084
#include "example_main.h"
#include <math.h>
#include <pil/pil_circuit.h>

#define V0 1.0
#define R_ 10.0
#define L_ 10e-3
#define C_ 10e-6
#define FA 100.0
#define FB 1000.0
#define TSWEEP 0.12

typedef struct {
  pil_celement *src;
  pil_celement *res;
  pil_celement *ind;
  pil_celement *cap;
} rlc_t;

static double drive_f(double t) { return FA + (FB - FA) * t / TSWEEP; }

static double chirp(double t, void *u) {
  (void)u;
  double beta = (FB - FA) / TSWEEP;

  return V0 * sin(2.0 * M_PI * (FA * t + 0.5 * beta * t * t));
}

static double sig_vs(const pil_world *w, void *u) {
  (void)w;

  return pil_celement_voltage(((rlc_t *)u)->src);
}

static double sig_vr(const pil_world *w, void *u) {
  (void)w;

  return R_ * pil_celement_current(((rlc_t *)u)->res);
}

static double sig_f(const pil_world *w, void *u) {
  (void)u;

  return drive_f(pil_world_time(w));
}

// Band-pass amplitude: R |I|, from i and di/dt = (v_s - R i - v_C) / L
static double sig_bp(const pil_world *w, void *u) {
  const rlc_t *r = u;
  double om = 2.0 * M_PI * drive_f(pil_world_time(w));
  double i = pil_celement_current(r->ind);
  double di =
      (pil_celement_voltage(r->src) - R_ * i - pil_celement_voltage(r->cap)) /
      L_;

  return R_ * sqrt(i * i + di * di / (om * om)) / V0;
}

// Low-pass amplitude: v_C and dv_C/dt = i / C
static double sig_lp(const pil_world *w, void *u) {
  const rlc_t *r = u;
  double om = 2.0 * M_PI * drive_f(pil_world_time(w));
  double vc = pil_celement_voltage(r->cap);
  double dv = pil_celement_current(r->ind) / C_;

  return sqrt(vc * vc + dv * dv / (om * om)) / V0;
}

int main(int argc, char **argv) {
  static rlc_t r;
  pil_world *w = pil_world_new();
  pil_circuit *c = pil_add_circuit(w);

  pil_node a = pil_circuit_node(c, 0.0, 1.6), g = pil_circuit_node(c, 0.0, 0.0),
           n1 = pil_circuit_node(c, 1.8, 1.6),
           n2 = pil_circuit_node(c, 3.6, 1.6),
           n3 = pil_circuit_node(c, 5.0, 1.6),
           br = pil_circuit_node(c, 5.0, 0.0);
  pil_circuit_ground(c, g);

  r.src = pil_vsource_fn(c, a, g, chirp, NULL, "v_s");
  r.res = pil_resistor(c, a, n1, R_, "R");
  r.ind = pil_inductor(c, n1, n2, L_, 0.0, "L");
  pil_wire(c, n2, n3);
  r.cap = pil_capacitor(c, n3, br, C_, 0.0, "C");
  pil_wire(c, br, g);
  pil_celement_set_label_at(r.src, -34, 0);
  pil_circuit_set_step(c, 2e-6);

  pil_plot *p = pil_add_plot(w, 0.0, 0.58, 0.46, 1.0, "t [s]", "[V]", TSWEEP);
  pil_plot_add_signal(p, "v_s", sig_vs, &r);
  pil_plot_add_signal(p, "v_R", sig_vr, &r);
  pil_plot *q = pil_add_plot(w, 0.52, 0.58, 0.99, 1.0, "f [Hz]",
                             "|v_{out}| / V_0", TSWEEP);
  pil_plot_add_signal(q, "", sig_f, NULL);
  pil_plot_add_signal(q, "band-pass (R)", sig_bp, &r);
  pil_plot_add_signal(q, "low-pass (C)", sig_lp, &r);
  pil_plot_set_xy(q, 0);
  pil_plots_record(w, 1e-4, 1e-4);

  pil_view view;
  pil_view_fit_region(&view, -0.8, -0.45, 5.6, 2.1, 1300, 760, 0.0, 0.0, 1.0,
                      0.57, 18);
  int rc = example_main(w, &view, "Frequency response of series RLC", argc,
                        argv, "bode.png");

  pil_world_free(w);

  return rc;
}
