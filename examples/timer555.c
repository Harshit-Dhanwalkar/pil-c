/*
 * timer555.c - 555 astable oscillator
 * C charges through R1 + R2 to 2/3 Vcc, discharges through R2 to 1/3 Vcc;
 * output is high while charging. R1 = 1 k, R2 = 2.2 k, C = 1 uF, Vcc = 9 V:
 * f = 1 / (ln2 (R1 + 2 R2) C) ~ 267 Hz, duty = (R1 + R2) / (R1 + 2 R2) ~ 59 %
 * Physics lives in pil_timer555; this file only plots it
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.002
#define PIL_EXAMPLE_IMAGE_TIME 0.02
#include "example_main.h"
#include <math.h>
#include <pil/pil_timer555.h>

#define R1 1000.0
#define R2 2200.0
#define C_ 1e-6
#define VCC 9.0
#define DUR 0.02

static pil_555 osc;
static double tnow = -1.0;

static void advance(const pil_world *w) {
  double t = pil_world_time(w);
  if (t < tnow || tnow < 0.0) {
    pil_555_init(&osc, R1, R2, C_, VCC);
    tnow = 0.0;
  }

  while (t - tnow > 0.0) {
    double h = t - tnow < 2e-5 ? t - tnow : 2e-5;
    pil_555_step(&osc, h);
    tnow += h;
  }
}

static double s_vc(const pil_world *w, void *u) {
  (void)u;
  advance(w);
  return osc.vc;
}
static double s_out(const pil_world *w, void *u) {
  (void)u;
  advance(w);
  return osc.out ? VCC : 0.0;
}
static double s_hi(const pil_world *w, void *u) {
  (void)w;
  (void)u;
  return 2.0 * VCC / 3.0;
}
static double s_lo(const pil_world *w, void *u) {
  (void)w;
  (void)u;
  return VCC / 3.0;
}

int main(int argc, char **argv) {
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);

  pil_plot *p =
      pil_add_plot(w, 0.08, 0.06, 0.98, 0.47, "t [s]", "v_C [V]", DUR);
  pil_plot_add_signal(p, "v_C", s_vc, NULL);
  pil_plot_add_signal(p, "2/3 V_{cc}", s_hi, NULL);
  pil_plot_add_signal(p, "1/3 V_{cc}", s_lo, NULL);
  pil_plot *q =
      pil_add_plot(w, 0.08, 0.56, 0.98, 0.97, "t [s]", "v_{out} [V]", DUR);
  pil_plot_add_signal(q, "v_{out}", s_out, NULL);
  pil_plots_record(w, 5e-5, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, 0, 0, 1, 1, 1000, 760, 0.0, 0.0, 1.0, 1.0, 0);
  int rc = example_main(w, &view, "555 astable", argc, argv, "timer555.png");

  pil_world_free(w);

  return rc;
}
