/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

#include "pil/pil_wave.h"
#include <math.h>

int pil_wave_received(const pil_wave_source *s, double xo, double t) {
  int n = 0;
  for (int i = 0; i * s->t0 <= t; i++) {
    double te = i * s->t0, xs = s->x0 + s->v * te;
    if (te + fabs(xo - xs) / s->c <= t) {
      n++;
    }
  }

  return n;
}

double pil_doppler_ratio(double c, double v, int ahead) {
  return ahead ? 1.0 / (1.0 - v / c) : 1.0 / (1.0 + v / c);
}
