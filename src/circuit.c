/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* circuit.c - lumped circuits by modified nodal analysis
 *
 * Unknowns: node voltages (all but ground) plus one branch current for every
 * voltage source and inductor. C and L use \theta-method companion model
 *   i_{n+1} = (C / (\theta h)) (v_{n+1} - v_n) - ((1-\theta)/\theta) i_n
 *   v_{n+1} = (L / (\theta h)) (i_{n+1} - i_n) - ((1-\theta)/\theta) v_n
 * (\theta = 1/2 is trapezoidal rule)
 * diode is linearised each Newton iteration with SPICE-style voltage limiting
 * Kirchhoff's current law is written as "sum of currents leaving a node = 0"
 */

#include "pil_circuit_priv.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define GMIN 1e-12
#define R_ON 1e-3
#define R_OFF 1e9

static void *grow(void *p, int n, int *cap, size_t sz) {
  if (n > *cap) {
    *cap = n * 2;
    p = realloc(p, (size_t)*cap * sz);
    if (!p) {
      abort();
    }
  }

  return p;
}

// Construction

pil_circuit *pil_add_circuit(pil_world *w) {
  pil_circuit *c = calloc(1, sizeof *c);
  if (!c) {
    abort();
  }

  c->world = w;
  c->max_h = 2e-5;
  c->theta = 0.5;
  c->first = 1;
  c->ground = -1;
  w->circuits =
      grow(w->circuits, w->ncircuits + 1, &w->ccap, sizeof *w->circuits);
  w->circuits[w->ncircuits++] = c;

  return c;
}

pil_node pil_circuit_node(pil_circuit *c, double x, double y) {
  c->nodes = grow(c->nodes, c->nn + 1, &c->ncap, sizeof *c->nodes);
  c->nodes[c->nn] = (pil_cnode){x, y, c->nn};

  return c->nn++;
}

void pil_circuit_ground(pil_circuit *c, pil_node n) { c->ground = n; }
void pil_circuit_set_step(pil_circuit *c, double max_h) { c->max_h = max_h; }
void pil_circuit_set_theta(pil_circuit *c, double theta) { c->theta = theta; }

int pil_cnode_find(const pil_circuit *c, int n) {
  while (c->nodes[n].parent != n) {
    n = c->nodes[n].parent;
  }

  return n;
}

void pil_wire(pil_circuit *c, pil_node a, pil_node b) {
  int ra = pil_cnode_find(c, a), rb = pil_cnode_find(c, b);
  if (ra != rb) {
    c->nodes[rb].parent = ra;
  }

  c->wires = grow(c->wires, c->nw + 1, &c->wcap, sizeof *c->wires);
  c->wires[c->nw++] = (pil_cwire){.a = a, .b = b};
}

static pil_celement *new_el(pil_circuit *c, pil_ce_kind k, pil_node a,
                            pil_node b, const char *label) {
  pil_celement *e = calloc(1, sizeof *e);
  if (!e) {
    abort();
  }

  e->kind = k;
  e->circuit = c;
  e->a = a;
  e->b = b;
  e->branch = -1;
  e->label_dy = 18;
  if (label) {
    snprintf(e->label, sizeof e->label, "%s", label);
  }

  c->els = grow(c->els, c->ne + 1, &c->ecap, sizeof *c->els);
  c->els[c->ne++] = e;

  return e;
}

pil_celement *pil_resistor(pil_circuit *c, pil_node a, pil_node b, double ohms,
                           const char *label) {
  pil_celement *e = new_el(c, CE_RESISTOR, a, b, label);
  e->value = ohms;

  return e;
}

pil_celement *pil_capacitor(pil_circuit *c, pil_node a, pil_node b,
                            double farads, double v0, const char *label) {
  pil_celement *e = new_el(c, CE_CAPACITOR, a, b, label);
  e->value = farads;
  e->init = v0;
  e->v = v0;

  return e;
}

pil_celement *pil_inductor(pil_circuit *c, pil_node a, pil_node b,
                           double henries, double i0, const char *label) {
  pil_celement *e = new_el(c, CE_INDUCTOR, a, b, label);
  e->value = henries;
  e->init = i0;
  e->i = i0;

  return e;
}

pil_celement *pil_vsource_dc(pil_circuit *c, pil_node plus, pil_node minus,
                             double volts, const char *label) {
  pil_celement *e = new_el(c, CE_VSOURCE, plus, minus, label);
  e->value = volts;

  return e;
}

pil_celement *pil_vsource_ac(pil_circuit *c, pil_node plus, pil_node minus,
                             double amp, double freq, double phase,
                             const char *label) {
  pil_celement *e = new_el(c, CE_VSOURCE, plus, minus, label);
  e->value = amp;
  e->freq = freq;
  e->phase = phase;
  e->ac = 1;

  return e;
}

pil_celement *pil_vsource_fn(pil_circuit *c, pil_node plus, pil_node minus,
                             pil_wave_fn fn, void *user, const char *label) {
  pil_celement *e = new_el(c, CE_VSOURCE, plus, minus, label);
  e->ac = 2;
  e->wave = fn;
  e->wave_user = user;

  return e;
}

pil_celement *pil_isource_dc(pil_circuit *c, pil_node from, pil_node to,
                             double amps, const char *label) {
  pil_celement *e = new_el(c, CE_ISOURCE, from, to, label);
  e->value = amps;

  return e;
}

pil_celement *pil_diode(pil_circuit *c, pil_node anode, pil_node cathode,
                        double is, double n_vt, const char *label) {
  pil_celement *e = new_el(c, CE_DIODE, anode, cathode, label);
  e->is = is;
  e->n_vt = n_vt;

  return e;
}

pil_celement *pil_switch(pil_circuit *c, pil_node a, pil_node b,
                         double t_change, int closed_before,
                         const char *label) {
  pil_celement *e = new_el(c, CE_SWITCH, a, b, label);
  e->t_change = t_change;
  e->closed_before = closed_before;

  return e;
}

static pil_celement *transducer(pil_circuit *c, pil_ce_kind k, pil_node a,
                                pil_node b, pil_body *m, pil_body *m2,
                                double kc, double r, double l,
                                const char *label) {
  pil_celement *e = new_el(c, k, a, b, label);
  e->mb = m;
  e->mb2 = m2;
  e->kc = kc;
  e->res = r;
  e->value = l;

  return e;
}

pil_celement *pil_motor(pil_circuit *c, pil_node a, pil_node b, pil_body *rotor,
                        pil_body *stator, double kt, double r, double l,
                        const char *label) {
  return transducer(c, CE_MOTOR, a, b, rotor, stator, kt, r, l, label);
}

pil_celement *pil_voice_coil(pil_circuit *c, pil_node a, pil_node b,
                             pil_body *mover, pil_body *base, pil_vec2 axis,
                             double kf, double r, double l, const char *label) {
  pil_celement *e = transducer(c, CE_VCOIL, a, b, mover, base, kf, r, l, label);
  double n = hypot(axis.x, axis.y);
  e->axis = n > 0 ? pil_v2(axis.x / n, axis.y / n) : pil_v2(1, 0);

  return e;
}

void pil_celement_set_label_at(pil_celement *e, double dx_px, double dy_px) {
  e->label_dx = dx_px;
  e->label_dy = dy_px;
}

void pil_celement_hide(pil_celement *e) { e->hidden = 1; }

// Results

double pil_celement_voltage(const pil_celement *e) { return e->v; }
double pil_celement_current(const pil_celement *e) { return e->i; }
double pil_celement_power(const pil_celement *e) { return e->v * e->i; }

double pil_circuit_node_voltage(const pil_circuit *c, pil_node n) {
  return c->vn ? c->vn[n] : 0.0;
}

int pil_circuit_newton_iterations(const pil_circuit *c) {
  return c->newton_iters;
}

double pil_circuit_node_current_sum(const pil_circuit *c, pil_node n) {
  int r = pil_cnode_find((pil_circuit *)c, n);
  double sum = 0;
  for (int k = 0; k < c->ne; k++) {
    const pil_celement *e = c->els[k];
    if (pil_cnode_find((pil_circuit *)c, e->a) == r) {
      sum += e->i;
    }

    if (pil_cnode_find((pil_circuit *)c, e->b) == r) {
      sum -= e->i;
    }
  }

  return sum;
}

double pil_circuit_energy(const pil_circuit *c) {
  double E = 0;
  for (int k = 0; k < c->ne; k++) {
    const pil_celement *e = c->els[k];
    if (e->kind == CE_CAPACITOR) {
      E += 0.5 * e->value * e->v * e->v;
    } else if (e->kind == CE_INDUCTOR || e->kind == CE_MOTOR ||
               e->kind == CE_VCOIL) {
      E += 0.5 * e->value * e->i * e->i;
    }
  }

  return E;
}

// Solver

static int solve_dense(double *A, double *b, int n) {
  for (int k = 0; k < n; k++) {
    int p = k;
    for (int i = k + 1; i < n; i++) {
      if (fabs(A[i * n + k]) > fabs(A[p * n + k])) {
        p = i;
      }
    }

    if (fabs(A[p * n + k]) < 1e-300) {
      return 0;
    }

    if (p != k) {
      for (int j = 0; j < n; j++) {
        double t = A[k * n + j];
        A[k * n + j] = A[p * n + j];
        A[p * n + j] = t;
      }

      double t = b[k];
      b[k] = b[p];
      b[p] = t;
    }

    for (int i = k + 1; i < n; i++) {
      double f = A[i * n + k] / A[k * n + k];
      if (f == 0.0) {
        continue;
      }

      for (int j = k; j < n; j++) {
        A[i * n + j] -= f * A[k * n + j];
      }

      b[i] -= f * b[k];
    }
  }

  for (int i = n - 1; i >= 0; i--) {
    double s = b[i];
    for (int j = i + 1; j < n; j++) {
      s -= A[i * n + j] * b[j];
    }

    b[i] = s / A[i * n + i];
  }

  return 1;
}

static double vsource_value(const pil_celement *e, double t) {
  if (e->ac == 1) {
    return e->value * sin(2.0 * M_PI * e->freq * t + e->phase);
  }

  if (e->ac == 2) {
    return e->wave(t, e->wave_user);
  }

  return e->value;
}

static int switch_closed(const pil_celement *e, double t) {
  return t >= e->t_change ? !e->closed_before : e->closed_before;
}

// SPICE pnjlim: keep exponential from running away
static double limit_diode(double vnew, double vold, double nvt, double vcrit) {
  if (vnew > vcrit && fabs(vnew - vold) > 2.0 * nvt) {
    if (vold > 0) {
      double arg = 1.0 + (vnew - vold) / nvt;

      return arg > 0 ? vold + nvt * log(arg) : vcrit;
    }

    return nvt * log(vnew / nvt);
  }

  return vnew;
}

// Assign wire currents from KCL. Element currents give external injection at
// each raw node; wires distribute flow through wire graph (assumed to be a
// forest; cycle edges get i = 0)
static void compute_wire_currents(pil_circuit *c) {
  int nn = c->nn;
  for (int k = 0; k < c->nw; k++) {
    c->wires[k].i = 0;
  }
  if (nn == 0 || c->nw == 0) {
    return;
  }

  double *ext = calloc(nn, sizeof *ext);
  double *acc = calloc(nn, sizeof *acc);
  int *pn = malloc(nn * sizeof *pn);
  int *pw = malloc(nn * sizeof *pw);
  int *bfs = malloc(nn * sizeof *bfs);
  char *seen = calloc(nn, 1);
  if (!ext || !acc || !pn || !pw || !bfs || !seen) {
    free(ext);
    free(acc);
    free(pn);
    free(pw);
    free(bfs);
    free(seen);
    return;
  }

  // External current leaving each raw node
  for (int k = 0; k < c->ne; k++) {
    const pil_celement *e = c->els[k];
    ext[e->a] += e->i;
    ext[e->b] -= e->i;
  }

  for (int start = 0; start < nn; start++) {
    if (seen[start]) {
      continue;
    }

    int qh = 0, qt = 0;
    bfs[qt++] = start;
    seen[start] = 1;
    pn[start] = -1;
    pw[start] = -1;

    while (qh < qt) {
      int u = bfs[qh++];
      for (int k = 0; k < c->nw; k++) {
        int wa = c->wires[k].a, wb = c->wires[k].b;
        int other = (wa == u) ? wb : ((wb == u) ? wa : -1);
        if (other < 0 || seen[other]) {
          continue;
        }
        seen[other] = 1;
        pn[other] = u;
        pw[other] = k;
        bfs[qt++] = other;
      }
    }

    // Reverse BFS: leaves first
    // F(u) = current leaving subtree u via its parent wire
    // KCL at u:  ext[u] + \sum F(children) - F(u) = 0
    for (int i = qt - 1; i >= 0; i--) {
      int u = bfs[i];
      double Fu = -ext[u] + acc[u];
      if (pw[u] >= 0) {
        int wk = pw[u], p = pn[u];
        int wa = c->wires[wk].a;
        c->wires[wk].i = (wa == u) ? Fu : -Fu;
        acc[p] += Fu;
      }
    }
  }

  free(ext);
  free(acc);
  free(pn);
  free(pw);
  free(bfs);
  free(seen);
}

static int step_once(pil_circuit *c, double h) {
  const double t1 = c->t + h;
  const double th = c->first ? 1.0 : c->theta;
  const double carry = (1.0 - th) / th;

  // node numbering: representatives only, ground excluded
  int *idx = malloc((size_t)c->nn * sizeof *idx);
  if (!idx) {
    abort();
  }

  int gr = c->ground >= 0 ? pil_cnode_find(c, c->ground) : -1, nv = 0;
  for (int n = 0; n < c->nn; n++) {
    idx[n] = -1;
  }

  for (int n = 0; n < c->nn; n++) {
    if (pil_cnode_find(c, n) == n && n != gr) {
      idx[n] = nv++;
    }
  }

  int nb = 0;
  for (int k = 0; k < c->ne; k++) {
    pil_celement *e = c->els[k];
    e->branch = -1;
    if (e->kind == CE_INDUCTOR || e->kind == CE_VSOURCE ||
        e->kind == CE_MOTOR || e->kind == CE_VCOIL) {
      e->branch = nv + nb++;
    }

    if (e->kind == CE_MOTOR) {
      e->wm = e->mb->omega - (e->mb2 ? e->mb2->omega : 0.0);
    } else if (e->kind == CE_VCOIL) {
      pil_vec2 vb = e->mb2 ? e->mb2->vel : pil_v2(0, 0);
      e->wm =
          (e->mb->vel.x - vb.x) * e->axis.x + (e->mb->vel.y - vb.y) * e->axis.y;
    }
  }

  int N = nv + nb;
  if (N > c->scap) {
    c->scap = N * 2;
    c->A = realloc(c->A, (size_t)c->scap * c->scap * sizeof *c->A);
    c->rhs = realloc(c->rhs, (size_t)c->scap * sizeof *c->rhs);
    c->x = realloc(c->x, (size_t)c->scap * sizeof *c->x);
    if (!c->A || !c->rhs || !c->x) {
      abort();
    }
  }

  double *A = c->A, *rhs = c->rhs;
  int has_diode = 0;
  for (int k = 0; k < c->ne; k++) {
    if (c->els[k]->kind == CE_DIODE) {
      has_diode = 1;
      c->els[k]->vd = c->els[k]->v; // start from last accepted voltage
    }
  }

#define IA(e) (idx[pil_cnode_find(c, (e)->a)])
#define IB(e) (idx[pil_cnode_find(c, (e)->b)])
#define AT(i, j) A[(size_t)(i) * N + (j)]

  int newton_ok = 0;
  int iter = 0;
  for (;; iter++) {
    memset(A, 0, (size_t)N * N * sizeof *A);
    memset(rhs, 0, (size_t)N * sizeof *rhs);
    for (int i = 0; i < nv; i++) {
      AT(i, i) += GMIN;
    }

    for (int k = 0; k < c->ne; k++) {
      pil_celement *e = c->els[k];
      int ia = IA(e), ib = IB(e);
      double G = 0, Ik = 0; // conductance and known current leaving a
      switch (e->kind) {
      case CE_RESISTOR: {
        G = 1.0 / e->value;
        break;
      }
      case CE_SWITCH: {
        G = 1.0 / (switch_closed(e, t1) ? R_ON : R_OFF);
        break;
      }
      case CE_CAPACITOR: {
        G = e->value / (th * h);
        Ik = -(G * e->v + carry * e->i);
        break;
      }
      case CE_DIODE: {
        double vd_c = e->vd;
        if (vd_c > 0.9) {
          vd_c = 0.9;
        }
        if (vd_c < -20.0) {
          vd_c = -20.0;
        }

        double ex = exp(fmin(vd_c / e->n_vt, 80.0));
        double id = e->is * (ex - 1.0);
        G = e->is * ex / e->n_vt + GMIN;
        if (G > 1e6) {
          G = 1e6;
        }
        Ik = id - G * vd_c;
        break;
      }
      case CE_ISOURCE: {
        Ik = e->value;
        break;
      }
      case CE_INDUCTOR: {
        double Rq = e->value / (th * h);
        if (ia >= 0) {
          AT(ia, e->branch) += 1;
          AT(e->branch, ia) += 1;
        }

        if (ib >= 0) {
          AT(ib, e->branch) -= 1;
          AT(e->branch, ib) -= 1;
        }

        AT(e->branch, e->branch) -= Rq;
        rhs[e->branch] = -Rq * e->i - carry * e->v;
        continue;
      }
      case CE_MOTOR:
      case CE_VCOIL: {
        // v = r i + l (i - i_n)/h + k w, backward Euler
        double Rq = e->res + e->value / h;
        if (ia >= 0) {
          AT(ia, e->branch) += 1;
          AT(e->branch, ia) += 1;
        }

        if (ib >= 0) {
          AT(ib, e->branch) -= 1;
          AT(e->branch, ib) -= 1;
        }

        AT(e->branch, e->branch) -= Rq;
        rhs[e->branch] = -(e->value / h) * e->i + e->kc * e->wm;
        continue;
      }
      case CE_VSOURCE: {
        if (ia >= 0) {
          AT(ia, e->branch) += 1;
          AT(e->branch, ia) += 1;
        }

        if (ib >= 0) {
          AT(ib, e->branch) -= 1;
          AT(e->branch, ib) -= 1;
        }

        rhs[e->branch] = vsource_value(e, t1);
        continue; // branch-row element: skip nodal stamping below
      }
      }

      // Nodal stamping for two-terminal conductance elements
      // (R, switch, C, diode, I-source).  L / V / motor / vcoil have already
      // `continue`d above; the cases above that `break` fall through to here.
      if (ia >= 0) {
        AT(ia, ia) += G;
        rhs[ia] -= Ik;
      }

      if (ib >= 0) {
        AT(ib, ib) += G;
        rhs[ib] += Ik;
      }

      if (ia >= 0 && ib >= 0) {
        AT(ia, ib) -= G;
        AT(ib, ia) -= G;
      }
    }

    if (!solve_dense(A, rhs, N)) {
      newton_ok = 0; // singular: keep previous state
      break;
    }

    if (!has_diode) {
      newton_ok = 1;
      break;
    }

    int converged = 1;
    for (int k = 0; k < c->ne; k++) {
      pil_celement *e = c->els[k];
      if (e->kind != CE_DIODE) {
        continue;
      }

      int ia = IA(e), ib = IB(e);
      double vsol = (ia >= 0 ? rhs[ia] : 0.0) - (ib >= 0 ? rhs[ib] : 0.0);
      double vcrit = e->n_vt * log(e->n_vt / (1.41421356 * e->is));
      double vn = limit_diode(vsol, e->vd, e->n_vt, vcrit);
      if (fabs(vsol - e->vd) > 1e-9 + 1e-6 * fabs(vsol)) {
        converged = 0;
      }

      e->vd = vn;
    }

    if (converged) {
      newton_ok = 1;
      break;
    }
    if (iter >= 99) {
      newton_ok = 0;
      break;
    }
  }

  c->newton_iters = iter + 1;

  // Accept: node voltages, then element states
  c->vn = realloc(c->vn, (size_t)c->nn * sizeof *c->vn);
  if (!c->vn) {
    abort();
  }

  for (int n = 0; n < c->nn; n++) {
    int i = idx[pil_cnode_find(c, n)];
    c->vn[n] = i >= 0 ? rhs[i] : 0.0;
  }

  for (int k = 0; k < c->ne; k++) {
    pil_celement *e = c->els[k];
    int ia = IA(e), ib = IB(e);
    double v = (ia >= 0 ? rhs[ia] : 0.0) - (ib >= 0 ? rhs[ib] : 0.0);
    double i = 0;
    switch (e->kind) {
    case CE_RESISTOR: {
      i = v / e->value;
      break;
    }
    case CE_SWITCH: {
      i = v / (switch_closed(e, t1) ? R_ON : R_OFF);
      break;
    }
    case CE_CAPACITOR: {
      double G = e->value / (th * h);
      i = G * (v - e->v) - carry * e->i;
      break;
    }
    case CE_DIODE: {
      double vd_c = e->vd;
      if (vd_c > 1.5) {
        vd_c = 1.5;
      }
      if (vd_c < -20.0) {
        vd_c = -20.0;
      }
      double ex = exp(fmin(vd_c / e->n_vt, 80.0));
      double g = e->is * ex / e->n_vt + GMIN;
      if (g > 1e6) {
        g = 1e6;
      }
      i = e->is * (ex - 1.0) + g * (v - vd_c);
      break;
    }
    case CE_ISOURCE: {
      i = e->value;
      break;
    }
    case CE_INDUCTOR:
    case CE_VSOURCE:
    case CE_MOTOR:
    case CE_VCOIL: {
      i = rhs[e->branch];
      break;
    }
    }

    e->v = v;
    e->i = i;
    // dots move in proportion to i / peak, so fastest branch always streams at
    // current_dot_speed whatever circuit's current scale
    if (fabs(i) > c->peak_i && isfinite(i)) {
      c->peak_i = fabs(i);
    }

    e->anim_phase += i / fmax(c->peak_i, 1e-12) * h;
  }

  compute_wire_currents(c);
  for (int k = 0; k < c->nw; k++) {
    c->wires[k].anim_phase += c->wires[k].i / fmax(c->peak_i, 1e-12) * h;
  }

#undef IA
#undef IB
#undef AT

  free(idx);

  c->t = t1;
  c->first = 0;

  return newton_ok ? 0 : 1;
}

static void circuit_step(pil_circuit *c, double dt) {
  if (c->nn == 0 || dt <= 0) {
    return;
  }

  // Very first step is backward Euler (C/L companion models need consistent
  // previous currents/voltages); make it so short that its damping is
  // negligible, then continue with chosen \theta
  if (c->first) {
    double h0 = fmin(c->max_h * 1e-3, 0.5 * dt);
    // if (step_once(c, h0) != 0) {
    //   // First step failed - skip it, keep going
    // } else {
    //   dt -= h0;
    // }
    step_once(c, h0); // ignore status on first step
    dt -= h0;
  }

  int n = (int)ceil(dt / c->max_h);
  if (n < 1) {
    n = 1;
  }

  // double h = dt / n;
  //
  // for (int s = 0; s < n; s++) {
  //   double hs = h;
  //   int tries = 0;
  //   while (step_once(c, hs) != 0) {
  //     // Save state, roll back and retry with half step
  //     pil_celement_state_save(c); // snapshot v, i, vd, vn, t, first
  //     hs *= 0.5;
  //     if (++tries >= 8) {
  //       pil_celement_state_restore(c); // give up, keep previous state
  //       break;
  //     }
  //   }
  // }

  for (int s = 0; s < n; s++) {
    step_once(c, dt / n);
  }
}

static void circuit_reset(pil_circuit *c) {
  c->t = 0;
  c->first = 1;
  c->newton_iters = 0;
  c->peak_i = 0;
  for (int k = 0; k < c->ne; k++) {
    pil_celement *e = c->els[k];
    e->v = e->kind == CE_CAPACITOR ? e->init : 0.0;
    e->i = e->kind == CE_INDUCTOR ? e->init : 0.0;
    e->anim_phase = 0;
    e->vd = 0;
  }

  for (int k = 0; k < c->nw; k++) {
    c->wires[k].i = 0;
    c->wires[k].anim_phase = 0;
  }

  free(c->vn);
  c->vn = NULL;
}

void pil_circuits_step(pil_world *w, double dt) {
  for (int i = 0; i < w->ncircuits; i++) {
    circuit_step(w->circuits[i], dt);
  }
}

void pil_circuits_reset(pil_world *w) {
  for (int i = 0; i < w->ncircuits; i++) {
    circuit_reset(w->circuits[i]);
  }
}

void pil_circuits_apply_forces(pil_world *w) {
  for (int n = 0; n < w->ncircuits; n++) {
    const pil_circuit *c = w->circuits[n];
    for (int k = 0; k < c->ne; k++) {
      const pil_celement *e = c->els[k];
      if (e->kind == CE_MOTOR) {
        pil_apply_torque(e->mb, e->kc * e->i);
        pil_apply_torque(e->mb2, -e->kc * e->i);
      } else if (e->kind == CE_VCOIL) {
        pil_vec2 f = pil_v2(e->kc * e->i * e->axis.x, e->kc * e->i * e->axis.y);
        pil_apply_force(e->mb, e->mb->pos, f);
        if (e->mb2) {
          pil_apply_force(e->mb2, e->mb2->pos, pil_v2(-f.x, -f.y));
        }
      }
    }
  }
}

void pil_circuits_free(pil_world *w) {
  for (int i = 0; i < w->ncircuits; i++) {
    pil_circuit *c = w->circuits[i];
    for (int k = 0; k < c->ne; k++) {
      free(c->els[k]);
    }

    free(c->els);
    free(c->nodes);
    free(c->wires);
    free(c->vn);
    free(c->A);
    free(c->rhs);
    free(c->x);
    free(c);
  }

  free(w->circuits);
  w->circuits = NULL;
  w->ncircuits = 0;
}
