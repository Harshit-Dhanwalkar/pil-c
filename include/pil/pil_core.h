/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/*
 * pil_core.h - physics core: 2D rigid bodies, force elements, constraints and a
 * DAE integrator
 *
 * Conventions: SI units, y axis up, angles in radians (counter-clockwise
 * positive), bodies are located by their centre of mass
 *
 * Use Cases:
 *
 *   pil_world *w = pil_world_new();
 *   pil_world_set_gravity(w, 0, -9.81);
 *   pil_body *m = pil_add_body(w, &(pil_body_desc){
 *        .shape = pil_rect(0.6, 0.3), .mass = 2.0, .pos = {0, 1}, .label = "m"
 * }); pil_add_spring(w, pil_at(0, 2), pil_on(m, 0, 0.15), 50.0,
 * PIL_REST_INITIAL); pil_world_advance(w, 1.0, 1e-3);
 */

#ifndef PIL_CORE_H
#define PIL_CORE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Forward declarations so core never needs <cairo.h>
 * (`struct _cairo` is type cairo calls cairo_t, so overlay callbacks can be
 * written with usual `cairo_t *`) */
struct _cairo;
typedef struct pil_view pil_view;

typedef struct pil_world pil_world;
typedef struct pil_element pil_element;
typedef struct pil_deco pil_deco;
typedef struct pil_terrain pil_terrain;

typedef struct {
  double x, y;
} pil_vec2;
static inline pil_vec2 pil_v2(double x, double y) {
  pil_vec2 v = {x, y};
  return v;
}

/* Shapes */
typedef enum {
  PIL_SHAPE_POINT = 0, /* point mass, drawn as a filled dot of radius r   */
  PIL_SHAPE_RECT,      /* w x h centred on shape origin               */
  PIL_SHAPE_DISC,      /* radius r                                        */
  PIL_SHAPE_POLY,      /* convex or concave polygon, local coordinates    */
  PIL_SHAPE_ARC_BAND   /* annular sector r_in..r_out, angles a0..a1       */
} pil_shape_kind;

typedef struct {
  pil_shape_kind kind;
  double w, h;                /* RECT                                      */
  double r;                   /* DISC / POINT                              */
  double r_in, r_out, a0, a1; /* ARC_BAND (angles in body frame)           */
  const pil_vec2 *pts;        /* POLY (copied by pil_add_body)             */
  int npts;
  pil_vec2 offset; /* shape origin relative to COM          */
} pil_shape;

pil_shape pil_point_mass(double draw_radius);
pil_shape pil_rect(double w, double h);
pil_shape pil_disc(double r);
pil_shape pil_poly(const pil_vec2 *pts, int n);
/* Annular sector (e.g. rolling "bowl" body).  If body's inertia is left at 0,
 * COM and inertia about it are computed automatically and shape.offset is set
 * so that body origin is COM */
pil_shape pil_arc_band(double r_in, double r_out, double a0, double a1);

/* Bodies */
enum {
  PIL_STATIC = 1 << 0,   /* infinite mass, does not move                 */
  PIL_HIDE_COM = 1 << 1, /* do not draw centre-of-mass symbol        */
  PIL_SPOKE = 1 << 2,    /* draw a radial tick on discs (shows rotation) */
  PIL_HIDDEN = 1 << 3,   /* do not draw body at all                  */
  PIL_NO_FILL = 1 << 4   /* outline only                                 */
};

typedef struct {
  pil_shape shape;
  double mass;    /* 0 -> 1 kg (ignored for pil_STATIC)            */
  double inertia; /* 0 -> computed from shape                  */
  pil_vec2 pos;   /* COM position                                  */
  double angle;
  pil_vec2 vel;
  double omega;
  unsigned flags;
  double fill;               /* grey level 0..1; 0 -> style default           */
  const char *label;         /* e.g. "M", "m_1"                               */
  double label_dx, label_dy; /* label offset from COM in pixels (dy up)*/
} pil_body_desc;

typedef struct pil_body {
  int id;
  pil_shape shape;
  double mass;
  double inertia;
  double inv_mass;
  double inv_inertia;
  unsigned flags;
  pil_vec2 pos;
  double angle;
  pil_vec2 vel;
  double omega;
  pil_vec2 acc;
  double alpha; /* from last solve                         */
  pil_vec2 force;
  double torque; /* accumulators, valid inside force fns   */
  double fill;
  char label[32];
  double label_dx;
  double label_dy;
  int dof; /* internal: -1 for static bodies           */
  pil_vec2 pos0;
  double angle0;
  pil_vec2 vel0;
  double omega0;
} pil_body;

pil_body *pil_add_body(pil_world *w, const pil_body_desc *d);
void pil_body_set_label(pil_body *b, const char *text, double dx_px,
                        double dy_px);
/* Set current *and* initial state (what pil_world_reset returns to). */
void pil_body_set_state(pil_body *b, pil_vec2 pos, double angle, pil_vec2 vel,
                        double omega);
pil_vec2 pil_body_point(const pil_body *b, pil_vec2 local); /* local -> world */
pil_vec2 pil_body_shape_origin(
    const pil_body
        *b); /* shape origin in body frame (arc-band circle centre) */
pil_vec2
pil_body_velocity_at(const pil_body *b,
                     pil_vec2 world_pt); /* velocity of material point */
void pil_apply_force(pil_body *b, pil_vec2 world_point,
                     pil_vec2 f); /* inside force callbacks */
void pil_apply_torque(pil_body *b, double torque);

/* Anchors: a point on a body, or a fixed point in world */
typedef struct {
  pil_body *body;
  pil_vec2 local;
} pil_anchor;
static inline pil_anchor pil_on(pil_body *b, double x, double y) {
  pil_anchor a = {b, {x, y}};
  return a;
}
static inline pil_anchor pil_at(double x, double y) {
  pil_anchor a = {0, {x, y}};
  return a;
}
pil_vec2 pil_anchor_world(pil_anchor a);
/* Anchor on `b` at point that is currently at world (x, y) - lets a figure be
 * authored in world coordinates whatever body's pose */
pil_anchor pil_on_at(pil_body *b, double x, double y);
/* Area centroid of a polygon: a pil_poly body's COM sits there, so .pos =
 * pil_poly_centroid(pts, n) places polygon exactly where its points say */
pil_vec2 pil_poly_centroid(const pil_vec2 *pts, int n);

/* World */
pil_world *pil_world_new(void);
void pil_world_free(pil_world *w);
void pil_world_set_gravity(pil_world *w, double gx, double gy);
/* Baumgarte constraint stabilisation (defaults 10, 10) */
void pil_world_set_stabilization(pil_world *w, double alpha, double beta);
/* Optional global viscous drag (1/s) on all dynamic bodies; default 0 */
void pil_world_set_drag(pil_world *w, double linear, double angular);
/* Relax to static equilibrium: integrates with time frozen and temporary
 * heavy damping, stores result as initial state */
void pil_world_settle(pil_world *w, double seconds);

void pil_world_step(pil_world *w, double dt); /* one RK4 step */
void pil_world_advance(pil_world *w, double duration, double max_dt);
double pil_world_time(const pil_world *w);
void pil_world_snapshot(pil_world *w); /* remember current state as initial */
void pil_world_reset(pil_world *w);

double pil_world_energy(const pil_world *w); /* kinetic + gravity + springs */
double
pil_world_constraint_error(const pil_world *w); /* max |\Phi| over all rows */

int pil_world_body_count(const pil_world *w);
pil_body *pil_world_body(pil_world *w, int index);

/* Mouse-style interaction: a soft spring from a body point to a target. */
pil_body *pil_world_pick(const pil_world *w, pil_vec2 world_pt,
                         pil_vec2 *local_out);
void pil_world_set_grab(pil_world *w, pil_body *b, pil_vec2 local,
                        pil_vec2 target, double k, double c);
void pil_world_clear_grab(pil_world *w);

/* Elements (everything can be "drop in") */
#define PIL_REST_INITIAL                                                       \
  (-1.0) /* spring rest length = distance at creation                          \
          */
#define PIL_LENGTH_INITIAL (-1.0)
enum { PIL_ROPE = 1 }; /* distance constraint: only resists stretching */

/* Force elements */
pil_element *pil_add_spring(pil_world *w, pil_anchor a, pil_anchor b, double k,
                            double rest);
pil_element *pil_add_damper(pil_world *w, pil_anchor a, pil_anchor b, double c);
/* Torsional spring+damper between two bodies' angles (NULL = world). */
pil_element *pil_add_torsion_spring(pil_world *w, pil_body *a, pil_body *b,
                                    double k, double rest_angle, double c);
pil_element *pil_add_force(pil_world *w, pil_anchor at,
                           pil_vec2 f); /* constant force */

typedef void (*pil_force_fn)(pil_world *w, double t, void *user);
/* Escape hatch: your callback may call pil_apply_force / pil_apply_torque. */
pil_element *pil_add_custom_force(pil_world *w, pil_force_fn fn, void *user);

/* Constraints */
pil_element *pil_add_revolute(pil_world *w, pil_anchor a,
                              pil_anchor b); /* pin joint  */
pil_element *pil_add_distance(pil_world *w, pil_anchor a, pil_anchor b,
                              double length, unsigned flags); /* rod / rope */
/* Rope over a fixed pulley: anchor A hangs from tangent point c1, anchor B
 * from tangent point c2 (world points; give same point twice for a
 * point-like pulley)
 * |A-c1| + |B-c2| = length and rope can only pull
 * For a pulley of radius R with vertical hanging segments use
 * c1 = centre - (R,0), c2 = centre + (R,0)
 * pulley itself is drawn by adding a static disc body */
pil_element *pil_add_pulley_rope(pil_world *w, pil_anchor a, pil_anchor b,
                                 pil_vec2 c1, pil_vec2 c2, double length);
/* Point `a` is confined to a line through `line_origin` (on body B or world)
 * with direction `dir` (in B's local frame)
 * Rotation stays free */
pil_element *pil_add_slot(pil_world *w, pil_anchor a, pil_anchor line_origin,
                          pil_vec2 dir);
pil_element *pil_add_fixed_angle(pil_world *w, pil_body *a, pil_body *b,
                                 double rel_angle);
/* Slot + fixed relative angle
 * Returns slot element */
pil_element *pil_add_prismatic(pil_world *w, pil_anchor a,
                               pil_anchor line_origin, pil_vec2 dir);
/* Prescribed constant angular velocity (wheels rolling on a moving road,
 * motors) */
pil_element *pil_add_rotation_drive(pil_world *w, pil_body *a, double omega);
/* Circle (centre at `center_local` on body, radius R) rolls without slipping on
 * a fixed straight line through `line_point` with direction `line_dir`
 * Body lies on left of direction of travel (line_dir = (1,0) -> body above
 * line)
 * Works for discs, cylinders and arc-band bodies (centre =
 * pil_body_shape_origin(b)) */
pil_element *pil_add_rolling(pil_world *w, pil_body *b, pil_vec2 center_local,
                             double radius, pil_vec2 line_point,
                             pil_vec2 line_dir);
/* A circle (centre anchor `a`, radius ra, fixed on rolling body) rolls without
 * slipping on a circle (centre anchor `b`, radius rb, fixed on track body or
 * world)
 * inside != 0: rolling circle runs on concave side (a ball in a bowl; pass
 * bowl's INNER radius and, for an arc band, shape origin as `b`)
 * inside == 0: it runs on convex side (a wheel on a cylinder)
 * Both bodies may move
 * Contact is bilateral: ball cannot leave surface */
pil_element *pil_add_circle_rolling(pil_world *w, pil_anchor a, double ra,
                                    pil_anchor b, double rb, int inside);
/* Linear kinematic link: qa - ratio * qb = const, where q is one coordinate
 * (x, y or angle) of a body (NULL = world, always 0)
 * Generic building lock for belts, gears, racks and rope-driven followers: with
 * it a pulley disc, a roller and a shuttle can be tied together with their
 * inertia */
typedef enum {
  PIL_COORD_X = 0,
  PIL_COORD_Y = 1,
  PIL_COORD_ANGLE = 2
} pil_coord;
pil_element *pil_add_link(pil_world *w, pil_body *a, pil_coord ca, pil_body *b,
                          pil_coord cb, double ratio);
/* Rigidly join two bodies at a point (pin + fixed relative angle)
 * Returns pin */
pil_element *pil_add_weld(pil_world *w, pil_anchor a, pil_anchor b);
/* Unilateral contact (penalty + regularised Coulomb friction) between corners /
 * rim of `body` and a line attached to `surface` (NULL = world), given by a
 * point and direction in surface frame
 * Body must be on left of direction
 * Use for blocks on wedges, cylinders on floors, dropped objects
 * Equal and opposite forces act on `surface` */
pil_element *pil_add_floor_contact(pil_world *w, pil_body *body,
                                   pil_anchor line_point, pil_vec2 line_dir,
                                   double k, double c, double mu);

/* Generic element properties */
void pil_element_set_label(pil_element *e, const char *text);
void pil_element_set_label_at(pil_element *e, const char *text, double dx_px,
                              double dy_px);
void pil_element_set_drawing(pil_element *e, int coils,
                             double width_px); /* springs/dampers */
void pil_element_hide(pil_element *e);
double pil_element_force(
    const pil_element *e); /* spring/damper force, constraint multiplier */

/* Terrain (road profile, ground): height as a function of s = x + v * t */
typedef double (*pil_height_fn)(double s, void *user);
pil_terrain *pil_add_terrain(pil_world *w, pil_height_fn fn, void *user,
                             double scroll_speed);
double pil_terrain_height(const pil_terrain *t, double x, double time);
/* Vertical tyre model: penalty spring k + damper c between disc and terrain */
pil_element *pil_add_terrain_contact(pil_world *w, pil_terrain *t,
                                     pil_body *wheel, double radius, double k,
                                     double c);

/* Annotations: purely visual */
/* Text supports  x_1  x_{12}  x^2  and UTF-8 (e.g. "\xcf\x86" for \phi) */
pil_deco *pil_add_label(pil_world *w, pil_anchor at, const char *text,
                        double dx_px, double dy_px);
pil_deco *pil_add_dimension(pil_world *w, pil_anchor a, pil_anchor b,
                            double offset_px, const char *text);
/* Arrow starting at `from`, pointing (dx,dy) pixels (dy up) */
pil_deco *pil_add_arrow(pil_world *w, pil_anchor from, double dx_px,
                        double dy_px, const char *text);
pil_deco *pil_add_line(pil_world *w, pil_anchor a, pil_anchor b, int dash_dot);
/* Leader line with arrowhead touching `at`; label sits at far end
 * (dx,dy pixels from `at`, dy up) - "M" pointer of textbook figures */
pil_deco *pil_add_leader(pil_world *w, pil_anchor at, double dx_px,
                         double dy_px, const char *text);
/* Angle arc at `center` from direction of p1 to direction of p2 (shorter way)
 * All three are anchors, so it follows moving bodies and scales with zoom */
pil_deco *pil_add_angle_between(pil_world *w, pil_anchor center, pil_anchor p1,
                                pil_anchor p2, double radius_px,
                                const char *text);
/* Hatched ground / wall: hatching is drawn on right of direction a->b */
pil_deco *pil_add_ground(pil_world *w, pil_vec2 a, pil_vec2 b);
/* Small "g" arrow in corner of textbook figures */
pil_deco *pil_add_gravity_indicator(pil_world *w, pil_vec2 at);
/* "g (x)": gravity pointing into page (plan-view figures) */
pil_deco *pil_add_gravity_into_page(pil_world *w, pil_vec2 at);
/* Arc from world angle `ref_angle` to live angle of `body` */
pil_deco *pil_add_angle_marker(pil_world *w, pil_anchor center,
                               double ref_angle, pil_body *body,
                               double radius_px, const char *text);
void pil_deco_hide(pil_deco *d);

/* Draw anything you like on top (called with view used for frame) */
typedef void (*pil_overlay_fn)(struct _cairo *cr, const pil_view *view,
                               void *user);
void pil_world_add_overlay(pil_world *w, pil_overlay_fn fn, void *user);

/* Plot panel: curves of any quantity of world against time, drawn in a
 * rectangle of figure (fractions of canvas, origin top-left)
 * Whole curve is computed up front by pil_plot_record(), then a cursor follows
 * animation mechanism
 * At most 6 signals per plot */
typedef struct pil_plot pil_plot;
typedef double (*pil_signal_fn)(const pil_world *w, void *user);
pil_plot *pil_add_plot(pil_world *w, double x0, double y0, double x1, double y1,
                       const char *xlabel, const char *ylabel, double t_span);
void pil_plot_add_signal(pil_plot *p, const char *label, pil_signal_fn fn,
                         void *user);
/* Run world from its current state for t_span (sampling every `sample_dt`,
 * integrating with steps <= max_dt), fill every plot, then restore state */
void pil_plots_record(pil_world *w, double sample_dt, double max_dt);
/* Read back recorded data (e.g. to export a CSV): sample count, time and value
 * of signal `sig` at sample `k` */
int pil_plot_samples(const pil_plot *p);
double pil_plot_time(const pil_plot *p, int k);
double pil_plot_value(const pil_plot *p, int sig, int k);

#ifdef __cplusplus
}
#endif
#endif /* pil_CORE_H */
