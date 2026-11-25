#include <math.h>
#include <pil/pil_biot.h>
#include <pil/pil_solenoid.h>

void pil_solenoid(double r, double L, int n_turns, double x, double z, int nseg,
                  double *bx, double *bz) {
  double sx = 0.0;
  double sz = 0.0;
  for (int k = 0; k < n_turns; k++) {
    double z0 = -0.5 * L + (k + 0.5) * L / n_turns;
    double ax;
    double az;
    pil_biot_loop(r, z0, x, z, nseg, &ax, &az);
    sx += ax;
    sz += az;
  }

  *bx = sx;
  *bz = sz;
}

double pil_solenoid_axis(double r, double L, double n, double z) {
  double a = z + 0.5 * L;
  double b = z - 0.5 * L;

  return 2.0 * M_PI * n * (a / sqrt(r * r + a * a) - b / sqrt(r * r + b * b));
}

double pil_solenoid_ideal(double n) { return 4.0 * M_PI * n; }
