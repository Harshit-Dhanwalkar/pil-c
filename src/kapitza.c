#include <math.h>
#include <pil/pil_kapitza.h>

void pil_kapitza_deriv(double t, const double *y, double *dydt, void *user) {
  const pil_kapitza *p = user;
  dydt[0] = y[1];
  dydt[1] = (p->g - p->a * p->w * p->w * cos(p->w * t)) / p->l * sin(y[0]) -
            p->gamma * y[1];
}

double pil_kapitza_pivot(const pil_kapitza *p, double t) {
  return p->a * cos(p->w * t);
}

double pil_kapitza_veff(const pil_kapitza *p, double th) {
  double s = sin(th);

  return -(p->g / p->l) * cos(th) +
         p->a * p->a * p->w * p->w / (4.0 * p->l * p->l) * s * s;
}

double pil_kapitza_critical_amp(double g, double l, double w) {
  return sqrt(2.0 * g * l) / w;
}

double pil_kapitza_slow_freq(const pil_kapitza *p) {
  double k = p->a * p->a * p->w * p->w / (2.0 * p->l * p->l) - p->g / p->l;

  return k > 0.0 ? sqrt(k) : 0.0;
}
