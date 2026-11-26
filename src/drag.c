#include <math.h>
#include <pil/pil_drag.h>
#include <pil/pil_ode.h>
#include <stddef.h>

void pil_drag_deriv(double t, const double *y, double *dydt, void *user) {
  (void)t;
  const pil_drag *p = user;
  double sp = hypot(y[2], y[3]);
  dydt[0] = y[2];
  dydt[1] = y[3];
  dydt[2] = -p->k * sp * y[2];
  dydt[3] = -p->g - p->k * sp * y[3];
}

double pil_drag_terminal(const pil_drag *p) {
  return p->k > 0.0 ? sqrt(p->g / p->k) : INFINITY;
}

int pil_drag_trace(const pil_drag *p, double h0, double v0, double theta,
                   double dt, double (*pts)[2], int max_pts, double *t_flight) {
  double y[4] = {0.0, h0, v0 * cos(theta), v0 * sin(theta)};
  double t = 0.0;
  int n = 0;
  pts[n][0] = y[0];
  pts[n++][1] = y[1];
  while (n < max_pts) {
    double px = y[0], py = y[1];
    pil_rk4_step(pil_drag_deriv, 4, t, y, dt, (void *)p);
    t += dt;
    if (y[1] < 0.0 && t > dt) {
      double f = py / (py - y[1]);
      pts[n][0] = px + f * (y[0] - px);
      pts[n++][1] = 0.0;
      t = t - dt + f * dt;
      break;
    }

    pts[n][0] = y[0];
    pts[n++][1] = y[1];
  }

  if (t_flight) {
    *t_flight = t;
  }

  return n;
}

double pil_drag_range(const pil_drag *p, double h0, double v0, double theta) {
  static double pts[4000][2];
  int n = pil_drag_trace(p, h0, v0, theta, 0.005, pts, 4000, NULL);

  return pts[n - 1][0];
}
