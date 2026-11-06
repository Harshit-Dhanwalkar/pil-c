/* Internal definitions shared by library sources */

#ifndef PIL_INTERNAL_H
#define PIL_INTERNAL_H

#include "pil/pil_core.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
  EL_SPRING,
  EL_DAMPER,
  EL_TORSION,
  EL_FORCE,
  EL_CUSTOM,
  EL_TERRAIN_CONTACT,
  EL_REVOLUTE,
  EL_DISTANCE,
  EL_SLOT,
  EL_FIXED_ANGLE,
  EL_ROT_DRIVE,
  EL_ROLLING,
  EL_CIRCLE_ROLL,
  EL_LINK,
  EL_PULLEY,
  EL_FLOOR,
  EL_STOP
} pil_el_kind;

struct pil_element {
  pil_el_kind kind;
  pil_world *world;
  pil_anchor a, b;
  double k, rest, c;   /* spring / damper / torsion parameters           */
  double length;       /* distance constraint                            */
  double phase;        /* fixed angle / drive phase / rolling offset     */
  double mu;           /* floor friction coefficient                     */
  int ca, cb;          /* link: coordinate of a / b (pil_coord)          */
  double radius2, ref; /* circle rolling: track radius, reference angle  */
  int sigma;           /* circle rolling: -1 inside track, +1 outside    */
  double omega;        /* drive rate                                     */
  pil_vec2 dir;        /* slot direction (unit, in B local frame)        */
  pil_vec2 p1, p2;     /* pulley rope: fixed tangent points              */
  pil_vec2 f;          /* constant force                                 */
  unsigned flags;
  pil_force_fn fn;
  void *user;
  pil_terrain *terrain;
  double radius;
  int rope_off; /* set by solver when a rope goes slack           */
  double last_force;
  /* drawing */
  int hidden;
  int coils;
  double width_px;
  char label[32];
  double label_dx, label_dy;
  int label_auto;
};

typedef enum {
  D_LABEL,
  D_DIM,
  D_ARROW,
  D_LINE,
  D_ANGLE,
  D_GROUND,
  D_GRAVITY,
  D_ANGLE_PTS,
  D_ENERGY_BAR
} pil_deco_kind;
struct pil_deco {
  pil_deco_kind kind;
  pil_anchor a, b, c;
  double dx;
  double dy;
  double offset;
  double radius;
  double ref;
  int dashed;
  int hidden;
  int reverse;
  pil_body *body;
  char text[48];
};

struct pil_terrain {
  pil_height_fn fn;
  void *user;
  double scroll;
};

typedef struct {
  pil_overlay_fn fn;
  void *user;
} pil_overlay;

/* \phi = position error, vel = d(\phi)/dt, ja/jb = Jacobian rows,
 * \gamma = RHS of J*a = \gamma (already includes Baumgarte terms) */
typedef struct {
  double phi, vel, ja[3], jb[3], gamma;
} pil_row;
typedef struct {
  pil_element *el;
  int idx;
  pil_row r;
  double lambda;
  double target;
} pil_rowrec;

typedef struct {
  char label[48];
  pil_signal_fn fn;
  void *user;
  double *v; /* samples, one per entry of plot's time axis */
} pil_signal;

#define PIL_PLOT_MAX_SIGNALS 6
struct pil_plot {
  double x0, y0, x1, y1; /* canvas fractions, origin top-left          */
  char xlabel[32], ylabel[32];
  double span; /* time axis runs 0 .. span                              */
  pil_signal sig[PIL_PLOT_MAX_SIGNALS];
  int nsig;
  double *t;
  int n, cap;
};

struct pil_world {
  pil_body **bodies;
  int nbodies, bcap, ndyn;
  pil_element **els;
  int nels, ecap;
  pil_deco **decos;
  int ndecos, dcap;
  pil_terrain **terrains;
  int nterr, tcap;
  pil_plot **plots;
  int nplots, plcap;
  struct pil_circuit **circuits;
  int ncircuits, ccap;
  pil_overlay *overlays;
  int noverlays;
  int ocap;
  pil_vec2 g;
  double t;
  double tcur;
  double alpha;
  double beta;
  double drag_lin;
  double drag_ang;
  double dissipated; /* kinetic energy lost in impacts (stops, rope snaps) */
  int frozen;
  int have_snapshot;
  struct {
    pil_body *b;
    pil_vec2 local, target;
    double k, c;
  } grab;
  /* scratch */
  double *K, *rhs;
  int kcap;
  pil_rowrec *rows;
  int rowcap;
  double *s0, *k1, *k2, *k3, *k4, *tmp;
  int scap;
};

/* vector helpers */
static inline pil_vec2 v_add(pil_vec2 a, pil_vec2 b) {
  return pil_v2(a.x + b.x, a.y + b.y);
}
static inline pil_vec2 v_sub(pil_vec2 a, pil_vec2 b) {
  return pil_v2(a.x - b.x, a.y - b.y);
}
static inline pil_vec2 v_mul(pil_vec2 a, double s) {
  return pil_v2(a.x * s, a.y * s);
}
static inline double v_dot(pil_vec2 a, pil_vec2 b) {
  return a.x * b.x + a.y * b.y;
}
static inline double v_cross(pil_vec2 a, pil_vec2 b) {
  return a.x * b.y - a.y * b.x;
}
static inline pil_vec2 v_perp(pil_vec2 a) { return pil_v2(-a.y, a.x); }
static inline pil_vec2 v_rot(pil_vec2 a, double th) {
  double c = cos(th), s = sin(th);
  return pil_v2(c * a.x - s * a.y, s * a.x + c * a.y);
}

/* anchor helpers: offset r from body COM in world axes, world point,
 * velocity */
pil_vec2 pil_anchor_r(pil_anchor a);
pil_vec2 pil_anchor_vel(pil_anchor a);

/* elements.c */
/* Fills up to 2 rows; returns row count. A slack rope returns 0 unless
 * force_active */
int pil_constraint_eval(const pil_element *e, double t, pil_row rows[2],
                        int force_active);
int pil_is_constraint(const pil_element *e);
/* Unilateral stop: fills row (phi <= 0 allowed) and returns its restitution */
double pil_stop_row(const pil_element *e, pil_row *r);
int pil_is_rope(const pil_element *e); /* one-sided (pull-only) constraint */
double pil_rope_length(const pil_element *e); /* current geometric length */
void pil_element_apply(pil_element *e, pil_world *w, double t);
double pil_element_energy(const pil_element *e, double t);
void pil_elements_rebase(
    pil_world *w); /* re-derive phases after a state reset */
static inline int pil_isdyn(const pil_body *b) { return b && b->dof >= 0; }

/* solver.c */
void pil_compute_accel(pil_world *w);
/* circuit.c (declared in pil_circuit_priv.h for circuit sources) */
void pil_circuits_apply_forces(pil_world *w);
void pil_circuits_step(pil_world *w, double dt);
void pil_circuits_reset(pil_world *w);
void pil_circuits_free(pil_world *w);
/* plot.c */
void pil_plots_free(pil_world *w);
/* inelastic velocity impulse when a rope snaps taut while still extending
 * (called after each step) */
void pil_project_rope_snaps(pil_world *w);

#endif
