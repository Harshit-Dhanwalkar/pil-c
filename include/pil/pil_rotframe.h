/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_rotframe.h - motion seen from a rotating frame: Coriolis effects */
#ifndef pil_ROTFRAME_H
#define pil_ROTFRAME_H

#ifdef __cplusplus
extern "C" {
#endif

/* Small-swing pendulum in a frame turning at rate omega_z about vertical
 * (omega_z = Omega sin(latitude) for a Foucault pendulum): with state
 * (x, y, vx, vy)
 *   x'' = -w0^2 x + 2 omega_z y',   y'' = -w0^2 y - 2 omega_z x'.
 * Swing plane turns at -omega_z. Use with pil_rk4_step (n = 4) */
typedef struct {
  double w0;      /* pendulum angular frequency sqrt(g / l) */
  double omega_z; /* frame rotation component along vertical */
} pil_foucault;
void pil_foucault_deriv(double t, const double *y, double *dydt, void *user);

/* Orientation of swing plane (mod pi) from second moments of samples x[i0..i1],
 * y[i0..i1] (use one full period for an exact result) */
double pil_plane_angle(const double *x, const double *y, int i0, int i1);

#ifdef __cplusplus
}
#endif

#endif
