/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_wave.h - wave sources and observers: retarded time, Doppler effect */
#ifndef pil_WAVE_H
#define pil_WAVE_H

#ifdef __cplusplus
extern "C" {
#endif

/* A source moving along x at constant speed v from x0, emitting a pulse every
 * period t0 starting at t = 0; waves travel at speed c
 * Number of pulses that an observer at xo has received by time t (pulse emitted
 * at te reaches observer at te + |xo - xs(te)| / c) */
typedef struct {
  double c, v, x0, t0;
} pil_wave_source;
int pil_wave_received(const pil_wave_source *s, double xo, double t);

/* Received frequency ratio f_obs / f_source for an observer ahead (ahead = 1)
 * or behind (ahead = 0) of a source moving at v: 1 / (1 -+ v/c) */
double pil_doppler_ratio(double c, double v, int ahead);

#ifdef __cplusplus
}
#endif

#endif
