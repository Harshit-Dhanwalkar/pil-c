/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/*
 * pil_circuit.h - lumped electric circuits that live in a pil_world
 *
 * A circuit is a set of nodes placed at points of figure and two-terminal
 * elements between them
 * Behind picture is modified nodal analysis with companion models for C and L
 * (\theta method, trapezoidal by default) and a Newton loop for diode
 * Circuit steps with world, so plots, recording, reset and animation window all
 * work as for mechanics
 *
 *   pil_circuit *c = pil_add_circuit(w);
 *   pil_node a = pil_circuit_node(c, 0, 1);
 *   pil_node b = pil_circuit_node(c, 2, 1);
 *   pil_node g = pil_circuit_node(c, 0, 0);
 *   pil_node h = pil_circuit_node(c, 2, 0);
 *   pil_circuit_ground(c, g);
 *   pil_vsource_dc(c, a, g, 5.0, "V");
 *   pil_resistor(c, a, b, 1e3, "R");
 *   pil_capacitor(c, b, h, 1e-6, 0.0, "C");
 *   pil_wire(c, g, h);
 *
 * Sign conventions: an element's voltage is v(a) - v(b); its current is current
 * flowing a -> b through it
 * A source drives current out of its positive terminal a
 */

#ifndef pil_CIRCUIT_H
#define pil_CIRCUIT_H

#include "pil_core.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct pil_circuit pil_circuit;
typedef struct pil_celement pil_celement;
typedef int pil_node;

/* a source waveform: value at time t */
typedef double (*pil_wave_fn)(double t, void *user);

pil_circuit *pil_add_circuit(pil_world *w);
/* A node drawn at world point (x, y) */
pil_node pil_circuit_node(pil_circuit *c, double x, double y);
/* 0 V reference (and a ground symbol is drawn there) */
void pil_circuit_ground(pil_circuit *c, pil_node n);
/* Internal step: largest integration step (default 20 microseconds) and \theta
 * (0.5 trapezoidal, default; 1 backward Euler for stiff switching) */
void pil_circuit_set_step(pil_circuit *c, double max_h);
void pil_circuit_set_theta(pil_circuit *c, double theta);

/* Ideal wire: two nodes are one electrical node; a line is drawn */
void pil_wire(pil_circuit *c, pil_node a, pil_node b);

/* Two-terminal elements. `label` is drawn beside symbol (sub/superscripts as in
 * all pil text); pass NULL for none */
pil_celement *pil_resistor(pil_circuit *c, pil_node a, pil_node b, double ohms,
                           const char *label);
/* v0: initial voltage v(a) - v(b) */
pil_celement *pil_capacitor(pil_circuit *c, pil_node a, pil_node b,
                            double farads, double v0, const char *label);
/* i0: initial current a -> b */
pil_celement *pil_inductor(pil_circuit *c, pil_node a, pil_node b,
                           double henries, double i0, const char *label);
pil_celement *pil_vsource_dc(pil_circuit *c, pil_node plus, pil_node minus,
                             double volts, const char *label);
/* amp * sin(2 pi freq t + phase) */
pil_celement *pil_vsource_ac(pil_circuit *c, pil_node plus, pil_node minus,
                             double amp, double freq, double phase,
                             const char *label);
pil_celement *pil_vsource_fn(pil_circuit *c, pil_node plus, pil_node minus,
                             pil_wave_fn fn, void *user, const char *label);
pil_celement *pil_isource_dc(pil_circuit *c, pil_node from, pil_node to,
                             double amps, const char *label);
/* Shockley diode, anode a -> cathode b
 * Typical silicon: is = 1e-12, n_vt = 0.026 (n = 1, 300 K) */
pil_celement *pil_diode(pil_circuit *c, pil_node anode, pil_node cathode,
                        double is, double n_vt, const char *label);
/* Switch: state `closed_before` until t_change, opposite afterwards */
pil_celement *pil_switch(pil_circuit *c, pil_node a, pil_node b,
                         double t_change, int closed_before, const char *label);

/* Electromechanical transducers: a winding (resistance r, inductance l) whose
 * current is tied to a body's motion
 *   motor:      torque kt * i on `rotor` (reaction on `stator`, NULL = world),
 *               back-EMF kt * (omega_rotor - omega_stator)
 *   voice coil: force kf * i along `axis` on `mover` (reaction on `base`),
 *               back-EMF kf * (v_mover - v_base) . axis
 * v(a) - v(b) = r i + l di/dt + back-EMF
 * Coupling is evaluated once per world step, so keep world step small against
 * mechanical dynamics */
pil_celement *pil_motor(pil_circuit *c, pil_node a, pil_node b, pil_body *rotor,
                        pil_body *stator, double kt, double r, double l,
                        const char *label);
pil_celement *pil_voice_coil(pil_circuit *c, pil_node a, pil_node b,
                             pil_body *mover, pil_body *base, pil_vec2 axis,
                             double kf, double r, double l, const char *label);

/* Placement of label: pixel offset from symbol centre (dy up) */
void pil_celement_set_label_at(pil_celement *e, double dx_px, double dy_px);
/* Do not draw this element (or its symbol; its wires stay) */
void pil_celement_hide(pil_celement *e);

/* Results (valid after each step; usable as plot signals). */
double pil_celement_voltage(const pil_celement *e); /* v(a) - v(b)       */
double pil_celement_current(const pil_celement *e); /* a -> b through it */
double pil_celement_power(const pil_celement *e);   /* absorbed v * i    */
double pil_circuit_node_voltage(const pil_circuit *c, pil_node n);
/* Kirchhoff's current law: element currents leaving electrical node `n` (all
 * nodes joined by wires) add up to zero
 * Returns that sum, a residual of a few 1e-12 A at most */
double pil_circuit_node_current_sum(const pil_circuit *c, pil_node n);
double pil_circuit_energy(const pil_circuit *c); /* stored in C and L       */
int pil_circuit_newton_iterations(const pil_circuit *c); /* of last step */

#ifdef __cplusplus
}
#endif

#endif
