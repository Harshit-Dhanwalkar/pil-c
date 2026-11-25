/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_solenoid.h - finite solenoid built from coaxial current loops
 * Units mu0 / 4 pi = 1 and unit current, as in pil_biot.h; axis is z */
#ifndef PIL_SOLENOID_H
#define PIL_SOLENOID_H

#ifdef __cplusplus
extern "C" {
#endif

/* Field at (x, 0, z) of n_turns loops of radius r spread uniformly over length
 * L centred at z = 0 (loop k at -L/2 + (k + 1/2) L / n_turns), each summed over
 * nseg segments
 * Returns Bx and Bz */
void pil_solenoid(double r, double L, int n_turns, double x, double z, int nseg,
                  double *bx, double *bz);

/* Exact on-axis field of continuous current sheet with n turns per unit length:
 * 2 pi n [ (z + L/2) / sqrt(r^2 + (z + L/2)^2) - (z - L/2) / sqrt(r^2 + (z -
 * L/2)^2) ] */
double pil_solenoid_axis(double r, double L, double n, double z);

/* Field inside an infinitely long solenoid: 4 pi n */
double pil_solenoid_ideal(double n);

#ifdef __cplusplus
}
#endif

#endif
