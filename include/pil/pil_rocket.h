/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_rocket.h - rocket equation with staging, vertical flight in uniform
 * gravity. During a burn  m dv/dt = mdot ve - m g  with propellant  mass
 * falling at rate mdot, so without gravity v = ve ln(m0 / m) (Tsiolkovsky)
 * When a stage's propellant is gone its dry mass is jettisoned and next stage
 * fires; after last stage rocket coasts */
#ifndef PIL_ROCKET_H
#define PIL_ROCKET_H

#ifdef __cplusplus
extern "C" {
#endif

#define PIL_ROCKET_MAX_STAGES 4

typedef struct {
  double m_dry;  /* Structure (dropped at burnout) */
  double m_fuel; /* Propellant */
  double mdot;   /* Propellant flow rate */
  double ve;     /* Exhaust speed */
} pil_stage;

typedef struct {
  double g;
  double payload;
  int n;
  pil_stage st[PIL_ROCKET_MAX_STAGES]; /* st[0] fires first */
} pil_rocket;

typedef struct {
  double t;
  double h;
  double v;
  int stage; /* Firing stage, n once all are spent */
  double mf; /* Propellant left in firing stage */
} pil_rocket_state;

void pil_rocket_init(pil_rocket_state *s);

/* Mass at this instant: payload + every stage not yet dropped + propellant */
double pil_rocket_mass(const pil_rocket *r, const pil_rocket_state *s);

/* True while an engine is firing */
int pil_rocket_burning(const pil_rocket *r, const pil_rocket_state *s);

/* Advance by dt (RK4); a burnout inside step is resolved exactly and remaining
 * time is spent in next stage or coasting */
void pil_rocket_step(const pil_rocket *r, pil_rocket_state *s, double dt);

/* Ideal velocity gain ve ln(m0 / m1) */
double pil_rocket_dv(double ve, double m0, double m1);

#ifdef __cplusplus
}
#endif

#endif
