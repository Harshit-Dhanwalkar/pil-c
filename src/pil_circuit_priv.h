/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_circuit_priv.h - circuit data structures shared by solver and renderer
 * (internal, not installed) */
#ifndef pil_CIRCUIT_PRIV_H
#define pil_CIRCUIT_PRIV_H

#include "pil/pil_circuit.h"
#include "pil_internal.h"

typedef enum {
  CE_RESISTOR,
  CE_CAPACITOR,
  CE_INDUCTOR,
  CE_VSOURCE,
  CE_ISOURCE,
  CE_DIODE,
  CE_SWITCH
} pil_ce_kind;

struct pil_celement {
  pil_ce_kind kind;
  pil_circuit *circuit;
  pil_node a, b;
  double value;       /* R, C, L, source amplitude                          */
  double freq, phase; /* ac source                                          */
  int ac;             /* vsource: 1 = sine, 0 = dc, 2 = user wave function  */
  pil_wave_fn wave;
  void *wave_user;
  double is, n_vt; /* diode                                                 */
  double t_change; /* switch                                                */
  int closed_before;
  double init; /* initial v (C) or i (L)                                    */
  double v, i; /* state after last step: v(a)-v(b), a->b current            */
  double vd;   /* diode linearisation point                                 */
  int branch;  /* MNA branch unknown for L and V (-1 if none)               */
  char label[32];
  double label_dx, label_dy;
  int hidden;
  double anim_phase; /* \int i dt — drives current-dot motion */
};

typedef struct {
  double x, y;
  int parent; /* union-find: wires merge nodes */
} pil_cnode;

typedef struct {
  pil_node a, b;
  double i;          /* current a -> b, from KCL pass                        */
  double anim_phase; /* \int i dt for dot animation                          */
} pil_cwire;

struct pil_circuit {
  pil_world *world;
  pil_cnode *nodes;
  int nn, ncap;
  pil_cwire *wires;
  int nw, wcap;
  pil_celement **els;
  int ne, ecap;
  pil_node ground;
  double max_h, theta, t;
  int first; /* first step after a reset is backward Euler */
  int newton_iters;
  double *vn; /* node voltages of last step, indexed by node id */
  /* scratch */
  double *A, *rhs, *x;
  int scap;
  double peak_i; /* running maximum |i|, for stable dot speed */
};

int pil_cnode_find(const pil_circuit *c, int n);

void pil_circuits_step(pil_world *w, double dt);
void pil_circuits_reset(pil_world *w);
void pil_circuits_free(pil_world *w);

#endif
