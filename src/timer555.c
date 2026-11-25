#include <math.h>
#include <pil/pil_timer555.h>

void pil_555_init(pil_555 *t, double r1, double r2, double c, double vcc) {
  t->r1 = r1;
  t->r2 = r2;
  t->c = c;
  t->vcc = vcc;
  t->vc = vcc / 3.0;
  t->out = 1;
}

void pil_555_step(pil_555 *t, double dt) {
  int guard = 0;
  while (dt > 0.0 && guard++ < 8) {
    double tau = (t->out ? t->r1 + t->r2 : t->r2) * t->c;
    double target = t->out ? t->vcc : 0.0;
    double thr = t->out ? 2.0 * t->vcc / 3.0 : t->vcc / 3.0;
    // Time to reach threshold from present voltage
    double tc = tau * log((target - t->vc) / (target - thr));
    if (tc > dt) {
      t->vc = target + (t->vc - target) * exp(-dt / tau);
      return;
    }

    if (tc < 0.0) {
      tc = 0.0;
    }
    t->vc = thr;
    t->out = !t->out;
    dt -= tc;
  }
}

double pil_555_high_time(double r1, double r2, double c) {
  return M_LN2 * (r1 + r2) * c;
}
double pil_555_low_time(double r1, double r2, double c) {
  (void)r1;
  return M_LN2 * r2 * c;
}
double pil_555_frequency(double r1, double r2, double c) {
  return 1.0 / (pil_555_high_time(r1, r2, c) + pil_555_low_time(r1, r2, c));
}
double pil_555_duty(double r1, double r2) {
  return (r1 + r2) / (r1 + 2.0 * r2);
}
