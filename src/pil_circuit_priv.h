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
  CE_SWITCH,
  CE_MOTOR,
  CE_VCOIL
} pil_ce_kind;

struct pil_celement {
  pil_ce_kind kind;
  pil_circuit *circuit;
  pil_node a;
  pil_node b;
  double value; /* R, C, L, source amplitude                          */
  double freq;  /* Frequency of AC source                             */
  double phase; /* Phase of AC source                                 */
  int ac;       /* vsource: 1 = sine, 0 = dc, 2 = user wave function  */
  pil_wave_fn wave;
  void *wave_user;
  double is;
  double n_vt;     /* Diode                                           */
  double t_change; /* Switch                                          */
  int closed_before;
  double init;   /* Initial v (C) or i (L)                            */
  pil_body *mb;  /* Transducer: driven body                           */
  pil_body *mb2; /* Transducer: reaction body                         */
  pil_vec2 axis; /* Voice coil direction (unit)                       */
  double kc;     /* Transducer constant                               */
  double res;    /* Transducer winding resistance                     */
  double wm;     /* Transducer: mechanical speed seen by winding      */
  double v;      /* State after last step: v(a)-v(b)                  */
  double i;      /* State after last step:  a->b current              */
  double vd;     /* Diode linearisation point                         */
  double vz;     /* zener breakdown voltage (0 = ordinary diode)      */
  struct pil_celement *mate; /* coupled inductor (transformer)        */
  double mutual; /* Mutual inductance with mate                       */
  int branch;    /* MNA branch unknown for L and V (-1 if none)       */
  char label[32];
  double label_dx;
  double label_dy;
  int hidden;
  double anim_phase; /* \int i dt — drives current-dot motion         */
};

typedef enum { DEV_OPAMP, DEV_NPN, DEV_PNP, DEV_NMOS, DEV_PMOS } pil_dev_kind;

struct pil_cdevice {
  pil_dev_kind kind;
  pil_circuit *circuit;
  pil_node n[3]; /* op-amp: + - out;  bjt: C B E; mosfet: D G S         */
  double vsat;   /* op-amp rail                                         */
  double beta;   /* bjt forward current gain; mosfet transconductance K */
  double vth;    /* mosfet threshold magnitude                          */
  double v1, v2; /* Linearisation: op-amp v+ - v-; bjt vbe, vbc (signed)*/
  double ti[3];  /* Current into each terminal after last step          */
  int branch;
  char label[32];
  double label_dx, label_dy;
};

typedef struct {
  double x;
  double y;
  int parent; /* Union-find: wires merge nodes */
} pil_cnode;

typedef struct {
  pil_node a;
  pil_node b;
  double i;          /* Current a -> b, from KCL pass                        */
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
  pil_cdevice **devs;
  int nd, dcap;
  pil_node ground;
  double max_h;
  double theta, t;
  int first; /* First step after a reset is backward Euler */
  int newton_iters;
  double *vn; /* Node voltages of last step, indexed by node id */
  /* Scratch */
  double *A;
  double *rhs;
  double *x;
  int scap;
  double peak_i; /* Running maximum |i|, for stable dot speed */
};

int pil_cnode_find(const pil_circuit *c, int n);

void pil_circuits_apply_forces(pil_world *w);
void pil_circuits_step(pil_world *w, double dt);
void pil_circuits_reset(pil_world *w);
void pil_circuits_free(pil_world *w);

#endif
