/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_emwave.h - plane electromagnetic waves in vacuum (eps0 = 1, so mu0 =
 * 1 / c^2; E and B in matching units) */
#ifndef pil_EMWAVE_H
#define pil_EMWAVE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Monochromatic plane wave travelling along +x with phase
 * phi = k x - omega t + phase, omega = c k:
 *   E = e0 (0, cos(phi), ellipticity sin(phi))
 *   B = (x_hat cross E) / c
 * ellipticity 0 = linear polarisation along y, +-1 = circular */
typedef struct {
  double c, e0, k, ellipticity, phase;
} pil_emwave;

/* Fields at position x and time t */
void pil_emwave_field(const pil_emwave *w, double x, double t, double E[3],
                      double B[3]);
/* Poynting vector S = E cross B / mu0 (energy flux) and energy density
 * u = eps0 E^2 / 2 + B^2 / (2 mu0) */
void pil_emwave_poynting(const pil_emwave *w, double x, double t, double S[3]);
double pil_emwave_energy_density(const pil_emwave *w, double x, double t);
/* Time-averaged intensity c eps0 <E^2> = c e0^2 (1 + ellipticity^2) / 2 */
double pil_emwave_intensity(const pil_emwave *w);
/* Radiation pressure on a perfect absorber, I / c and on a perfect
 * reflector, 2 I / c */
double pil_emwave_pressure(const pil_emwave *w, int reflector);

#ifdef __cplusplus
}
#endif

#endif
