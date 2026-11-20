/*
 * freebody.c - torque-free asymmetric rigid body, I_1 < I_2 < I_3
 * Three bars cross at centre of mass
 * With no torque L = I w stays fixed in space (dash-dot axis) and so does
 * E_{rot} = 1/2 w . L, while angular velocity w wanders in body: spinning about
 * intermediate axis e_2 is unstable, body tumbles and w_2 changes sign
 * Centre of mass drifts along L at constant speed, so tip of long bar draws a
 * wake Euler's equations are solved with pil_rigid3 (z-x-z angles, L along
 * space z); drawing is in an overlay
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 7.0
#include "example_main.h"
#include <math.h>
#include <pil/pil_lagrange.h>

#define NTR 600
#define TR_DT 0.02

typedef struct {
  pil_world *w;
  pil_lsys *s;
  pil_rigid3 *p;
  pil_proj3 pr;
  double a[3];        // Half-lengths of bars along e_1, e_2, e_3
  double drift;       // Speed of CM along L
  double x0;          // CM position along L at t = 0
  pil_vec2 wake[NTR]; // Tip of bar 1, relative to CM, at samples
  double wakeP[NTR][3];
  int nwake;
} scene;

// Scene axes: proj x = space z (L axis, horizontal), proj y = space x, proj
// z = space y (up on screen)
static void to_proj(const double s[3], double out[3]) {
  out[0] = s[2];
  out[1] = s[0];
  out[2] = s[1];
}

static void body_to_space(const double R[3][3], const double b[3],
                          double out[3]) {
  for (int i = 0; i < 3; i++) {
    out[i] = R[i][0] * b[0] + R[i][1] * b[1] + R[i][2] * b[2];
  }
}

static pil_vec2 scr(const scene *sc, const pil_view *v, const double P[3]) {
  return pil_view_to_screen(v, pil_proj3_point(&sc->pr, P[0], P[1], P[2]));
}

static double cm_x(const scene *sc, double t) { return sc->x0 + sc->drift * t; }

static void precompute(scene *sc, const double w0[3]) {
  pil_world *w2 = pil_world_new();
  pil_lsys *s2 = pil_add_lagrangian(w2, 3, pil_rigid3_T, NULL, sc->p);
  double q0[3];
  double qd0[3];
  pil_rigid3_state(sc->p, w0, 0.0, q0, qd0);
  pil_lsys_set_state(s2, q0, qd0);
  for (int k = 0; k < NTR; k++) {
    double R[3][3];
    pil_rigid3_axes(pil_lsys_q(s2), R);
    double tip_b[3] = {sc->a[0], 0, 0};
    double tip_s[3];
    body_to_space(R, tip_b, tip_s);
    to_proj(tip_s, sc->wakeP[k]);
    sc->wakeP[k][0] += cm_x(sc, k * TR_DT);
    for (int i = 0; i < 20; i++) {
      pil_world_step(w2, TR_DT / 20.0);
    }
  }

  sc->nwake = NTR;
  pil_world_free(w2);
}

typedef struct {
  double v[4][3];
} face;

// Overlay to loop animation in --live mode
static void loop_overlay(cairo_t *cr, const pil_view *v, void *user) {
  (void)cr;
  (void)v;
  pil_world *w = user;
  // 600 points * 0.02 s = 12 seconds of simulated time
  if (pil_world_time(w) >= 12.0) {
    pil_world_reset(w);
  }
}

static void draw_overlay(cairo_t *cr, const pil_view *v, void *user) {
  scene *sc = user;
  pil_style st = pil_style_textbook();
  double t = pil_world_time(sc->w);
  const double *q = pil_lsys_q(sc->s);
  const double *qd = pil_lsys_qd(sc->s);
  double R[3][3];
  double wb[3];
  double ws[3];
  double wp[3];
  pil_rigid3_axes(q, R);
  pil_rigid3_omega(q, qd, wb);
  body_to_space(R, wb, ws);
  to_proj(ws, wp);
  double X = cm_x(sc, t);
  double Lm = 0;
  {
    double Lb[3] = {sc->p->I[0] * wb[0], sc->p->I[1] * wb[1],
                    sc->p->I[2] * wb[2]};
    Lm = sqrt(Lb[0] * Lb[0] + Lb[1] * Lb[1] + Lb[2] * Lb[2]);
  }

  cairo_save(cr);
  cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);

  // Axis of L through the CM, dash-dot
  {
    double p0[3] = {sc->x0 - 0.5, 0, 0};
    double p1[3] = {X + 2.0, 0, 0};
    pil_vec2 a = scr(sc, v, p0);
    pil_vec2 b = scr(sc, v, p1);
    double dd[] = {10, 3, 2, 3};
    cairo_set_line_width(cr, st.thin_w);
    cairo_set_dash(cr, dd, 4, 0);
    cairo_move_to(cr, a.x, a.y);
    cairo_line_to(cr, b.x, b.y);
    cairo_stroke(cr);
    cairo_set_dash(cr, NULL, 0, 0);
  }

  // Wake of bar tip, older part lighter
  int n = (int)(t / TR_DT);
  if (n > sc->nwake - 1) {
    n = sc->nwake - 1;
  }

  cairo_set_line_width(cr, 1.8);
  for (int k = 1; k <= n; k++) {
    double g = 0.1 + 0.65 * (1.0 - (double)k / (double)(n + 1));
    cairo_set_source_rgb(cr, g, g, g);
    pil_vec2 a = scr(sc, v, sc->wakeP[k - 1]), b = scr(sc, v, sc->wakeP[k]);
    cairo_move_to(cr, a.x, a.y);
    cairo_line_to(cr, b.x, b.y);
    cairo_stroke(cr);
  }

  // Bars as boxes: six faces each, painted far to near
  struct {
    face f;
    double depth;
    double fill;
  } faces[18];
  int nf = 0;
  const double th = 0.075;
  for (int bar = 0; bar < 3; bar++) {
    double h[3] = {th, th, th};
    h[bar] = sc->a[bar];
    for (int ax = 0; ax < 3; ax++) {
      for (int sgn = -1; sgn <= 1; sgn += 2) {
        int u = (ax + 1) % 3;
        int w = (ax + 2) % 3;
        double cs[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
        for (int c = 0; c < 4; c++) {
          double b[3];
          b[ax] = sgn * h[ax];
          b[u] = cs[c][0] * h[u];
          b[w] = cs[c][1] * h[w];
          double s[3];
          double P[3];
          body_to_space(R, b, s);
          to_proj(s, P);
          P[0] += X;
          for (int k = 0; k < 3; k++) {
            faces[nf].f.v[c][k] = P[k];
          }
        }

        double mid[3] = {0, 0, 0};
        for (int c = 0; c < 4; c++) {
          for (int k = 0; k < 3; k++) {
            mid[k] += 0.25 * faces[nf].f.v[c][k];
          }
        }

        faces[nf].depth = pil_proj3_depth(&sc->pr, mid[0], mid[1], mid[2]);
        faces[nf].fill = 0.72 + 0.1 * (bar == 1);
        nf++;
      }
    }
  }

  for (int i = 1; i < nf; i++) {
    for (int j = i; j > 0 && faces[j - 1].depth > faces[j].depth; j--) {
      __typeof__(faces[0]) tmp = faces[j];
      faces[j] = faces[j - 1];
      faces[j - 1] = tmp;
    }
  }

  cairo_set_line_width(cr, st.thin_w + 0.3);
  for (int i = 0; i < nf; i++) {
    for (int c = 0; c < 4; c++) {
      pil_vec2 p = scr(sc, v, faces[i].f.v[c]);
      if (c == 0) {
        cairo_move_to(cr, p.x, p.y);
      } else {
        cairo_line_to(cr, p.x, p.y);
      }
    }

    cairo_close_path(cr);
    cairo_set_source_rgb(cr, faces[i].fill, faces[i].fill, faces[i].fill);
    cairo_fill_preserve(cr);
    cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);
    cairo_stroke(cr);
  }

  // Body axes e_i, omega and L as arrows from CM
  pil_vec2 C = scr(sc, v, (double[3]){X, 0, 0});
  const char *names[3] = {"e_1", "e_2", "e_3"};
  const double elen[3] = {1.45, 1.3, 1.2};
  for (int i = 0; i < 3; i++) {
    double b[3] = {0, 0, 0}, s[3], P[3];
    b[i] = elen[i];
    body_to_space(R, b, s);
    to_proj(s, P);
    P[0] += X;
    pil_vec2 T = scr(sc, v, P);
    cairo_set_line_width(cr, st.line_w);
    pil_draw_arrow(cr, C, T, 10, &st);
    pil_draw_text(cr, (pil_vec2){T.x + (T.x > C.x ? 6 : -6), T.y - 8}, names[i],
                  T.x > C.x ? 0.0 : 1.0, &st);
  }

  {
    const double k = 0.45;
    double P[3] = {X + k * wp[0], k * wp[1], k * wp[2]};
    pil_vec2 T = scr(sc, v, P);
    cairo_set_line_width(cr, st.line_w + 1.5);
    pil_draw_arrow(cr, C, T, 12, &st);
    pil_draw_text(cr, (pil_vec2){T.x + 6, T.y + 18}, "\xcf\x89", 0.0, &st);
    double Lp[3] = {X + 0.62 * Lm / (Lm > 0 ? Lm : 1) * 2.2, 0, 0};
    pil_vec2 Lt = scr(sc, v, Lp);
    cairo_set_line_width(cr, st.thin_w);
    pil_draw_arrow(cr, C, Lt, 8, &st);
    pil_draw_text(cr, (pil_vec2){Lt.x + 4, Lt.y - 14}, "L", 0.0, &st);
  }

  // CM marker: quartered disc
  cairo_set_line_width(cr, st.thin_w);
  cairo_new_path(cr);
  cairo_arc(cr, C.x, C.y, 6, 0, 2 * M_PI);
  cairo_set_source_rgb(cr, 1, 1, 1);
  cairo_fill_preserve(cr);
  cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);
  cairo_stroke(cr);
  cairo_move_to(cr, C.x, C.y);
  cairo_arc(cr, C.x, C.y, 6, 0, 0.5 * M_PI);
  cairo_close_path(cr);
  cairo_fill(cr);
  cairo_move_to(cr, C.x, C.y);
  cairo_arc(cr, C.x, C.y, 6, M_PI, 1.5 * M_PI);
  cairo_close_path(cr);
  cairo_fill(cr);
  pil_draw_text(cr, (pil_vec2){C.x + 8, C.y - 12}, "CM", 0.0, &st);

  // Equations
  pil_draw_text(cr, (pil_vec2){v->width * 0.30, 52},
                "\xcf\x84 = 0,      I_1 < I_2 < I_3", 0.0, &st);
  pil_draw_text(
      cr, (pil_vec2){v->width * 0.12, 98},
      "L = I \xcf\x89 = const,      E_{rot} = 1/2 \xcf\x89 \xc2\xb7 L = "
      "const",
      0.0, &st);
  cairo_restore(cr);
}

static double sig_w(const pil_world *w, void *user, int i) {
  (void)w;
  scene *sc = user;
  double wb[3];
  pil_rigid3_omega(pil_lsys_q(sc->s), pil_lsys_qd(sc->s), wb);

  return wb[i];
}

static double sig_w1(const pil_world *w, void *u) { return sig_w(w, u, 0); }
static double sig_w2(const pil_world *w, void *u) { return sig_w(w, u, 1); }
static double sig_w3(const pil_world *w, void *u) { return sig_w(w, u, 2); }

int main(int argc, char **argv) {
  static pil_rigid3 p = {{1.0, 1.36, 1.64}};
  static scene sc;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);

  sc.w = w;
  sc.p = &p;
  sc.a[0] = 0.95;
  sc.a[1] = 0.76;
  sc.a[2] = 0.57;
  sc.drift = 0.45;
  sc.x0 = -1.2;
  sc.pr = pil_proj3_make(0, 0, 0.55, 0.30, 1.0);

  double w0[3] = {0.15, 3.0, 0.12};
  double q0[3], qd0[3];
  pil_rigid3_state(&p, w0, 0.4, q0, qd0);
  sc.s = pil_add_lagrangian(w, 3, pil_rigid3_T, NULL, &p);
  pil_lsys_set_state(sc.s, q0, qd0);
  precompute(&sc, w0);
  
  pil_world_add_overlay(w, loop_overlay, w);
  pil_world_add_overlay(w, draw_overlay, &sc);

  pil_plot *pl = pil_add_plot(w, 0.08, 0.76, 0.98, 0.97, "t [s]",
                              "\xcf\x89_i [rad/s]", 12.0);
  pil_plot_add_signal(pl, "\xcf\x89_1", sig_w1, &sc);
  pil_plot_add_signal(pl, "\xcf\x89_2", sig_w2, &sc);
  pil_plot_add_signal(pl, "\xcf\x89_3", sig_w3, &sc);
  pil_plots_record(w, 0.02, 5e-4);

  pil_view view;
  pil_view_fit_region(&view, -2.0, -1.7, 6.0, 1.7, 1500, 1000, 0.0, 0.12, 1.0,
                      0.74, 12);
  int rc = example_main(w, &view, "Torque-free rigid body", argc, argv,
                        "freebody.png");

  pil_world_free(w);

  return rc;
}
