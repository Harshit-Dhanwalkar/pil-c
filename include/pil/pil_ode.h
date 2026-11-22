/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_ode.h - fixed-step integrator for small ODE systems */
#ifndef pil_ODE_H
#define pil_ODE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Right-hand side of y' = f(t, y) for n unknowns */
typedef void (*pil_ode_fn)(double t, const double *y, double *dydt, void *user);

/* One classical fourth-order Runge-Kutta step of size h, y updated in place
 * (n <= 16) */
void pil_rk4_step(pil_ode_fn f, int n, double t, double *y, double h,
                  void *user);

/* Angle wrapped into [-pi, pi) */
double pil_wrap_pi(double a);

#ifdef __cplusplus
}
#endif

#endif
