/*
 * linkage.c - parallelogram suspension carrying a rolling disc, a slider and a
 * double pendulum
 * Bar 1 hangs from two equal massless rods of length l (angle phi) and so only
 * translates; a torsion spring k_phi acts on left rod. On bar a disc 2
 * rolls against spring k_1 and dashpot c_1; a block 3 slides in a slot of
 * disc against spring k_2; block 6 rests on bar with friction
 * Under bar a double pendulum (rods 7 and 8, ball 9) hangs from S_1 Five
 * coordinates: phi, disc position, slider position, theta_7, theta_8 Everything
 * is public-API
 */

#include "example_main.h"
#include <math.h>

typedef struct {
  pil_body *bar;
  pil_vec2 pin; // left pin of bar, in bar-frame offset from COM
  pil_vec2 pivot;
  double l, kphi;
} lk_ref;

typedef struct {
  pil_body *bar;
  pil_body *disc;
  pil_body *blk;
  pil_body *r7;
  pil_body *r8;
  pil_vec2 pin; // left pin of bar relative to its COM
  pil_vec2 pivot;
  pil_vec2 slot; // slot direction in disc frame
  double wall;   // wall x offset on bar
} lk_sig;

static pil_vec2 left_pin(const lk_sig *s) {
  return pil_body_point(s->bar, s->pin);
}

static double sig_phi(const pil_world *w, void *user) {
  const lk_sig *s = user;
  pil_vec2 p = left_pin(s);
  (void)w;

  return atan2(p.x - s->pivot.x, s->pivot.y - p.y);
}

static double sig_disc(const pil_world *w, void *user) {
  const lk_sig *s = user;
  (void)w;

  return s->disc->pos.x - left_pin(s).x;
}

static double sig_slot(const pil_world *w, void *user) {
  const lk_sig *s = user;
  pil_vec2 u = {
      s->slot.x * cos(s->disc->angle) - s->slot.y * sin(s->disc->angle),
      s->slot.x * sin(s->disc->angle) + s->slot.y * cos(s->disc->angle)};
  (void)w;

  return (s->blk->pos.x - s->disc->pos.x) * u.x +
         (s->blk->pos.y - s->disc->pos.y) * u.y;
}

static double sig_th7(const pil_world *w, void *user) {
  (void)w;

  return ((const lk_sig *)user)->r7->angle + M_PI / 2;
}

static double sig_th8(const pil_world *w, void *user) {
  (void)w;

  return ((const lk_sig *)user)->r8->angle + M_PI / 2;
}

// Torsion spring k_phi on left rod: it acts on rod angle, and rods are
// massless, so it is applied as equivalent force on bar's pin (virtual work Q
// dphi = F . l t dphi)
static void torsion_force(pil_world *w, double t, void *user) {
  const lk_ref *r = user;
  pil_vec2 p = pil_body_point(r->bar, r->pin);
  double phi = atan2(p.x - r->pivot.x, r->pivot.y - p.y);
  double F = -r->kphi * phi / r->l;
  (void)w;
  (void)t;

  pil_apply_force(r->bar, p, (pil_vec2){F * cos(phi), F * sin(phi)});
}

typedef struct {
  const lk_sig *s;
  double kphi_r; // radius of spiral, in px
  double a, l;
} lk_draw;

static void pin_symbol(cairo_t *cr, const pil_view *v, pil_vec2 p,
                       const pil_style *st) {
  pil_draw_hinge(cr, pil_view_to_screen(v, p), st);
}

static void lk_overlay(cairo_t *cr, const pil_view *v, void *user) {
  const lk_draw *d = user;
  const lk_sig *s = d->s;
  pil_style st = pil_style_textbook();
  pil_vec2 L = left_pin(s);
  pil_vec2 P1 = s->pivot;
  pil_vec2 P2 = {s->pivot.x + d->a, s->pivot.y};
  pil_vec2 Rr = {L.x + d->a, L.y};

  cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);
  pin_symbol(cr, v, P1, &st);
  pin_symbol(cr, v, P2, &st);
  pin_symbol(cr, v, L, &st);
  pin_symbol(cr, v, Rr, &st);

  // Torsion spring k_phi at left pivot: small spiral
  {
    pil_vec2 c = pil_view_to_screen(v, P1);
    cairo_set_line_width(cr, st.thin_w);
    for (int i = 0; i <= 90; i++) {
      double t = 4.5 * M_PI * i / 90.0;
      double r = 6 + 1.9 * t;
      double x = c.x + r * cos(t);
      double y = c.y + r * sin(t);
      if (i == 0) {
        cairo_move_to(cr, x, y);
      } else {
        cairo_line_to(cr, x, y);
      }
    }
    cairo_stroke(cr);
    pil_draw_text(cr, (pil_vec2){c.x + 30, c.y + 26}, "k_\xcf\x86", 0, &st);
  }

  // Dashpot c_1 beside spring k_1: wall to disc centre, drawn offset
  {
    pil_vec2 a = {L.x + s->wall, s->disc->pos.y};
    pil_vec2 b = s->disc->pos;
    pil_vec2 sa = pil_view_to_screen(v, a);
    pil_vec2 sb = pil_view_to_screen(v, b);
    sa.y += 30;
    sb.y += 30;
    cairo_set_line_width(cr, st.thin_w);
    pil_draw_dashpot(cr, sa, sb, 150.0, 7, &st);
    cairo_move_to(cr, sa.x, sa.y);
    cairo_line_to(cr, sa.x, sa.y - 30);
    cairo_move_to(cr, sb.x, sb.y);
    cairo_line_to(cr, sb.x, sb.y - 30);
    cairo_stroke(cr);
    pil_draw_text(cr, (pil_vec2){(sa.x + sb.x) / 2, sa.y + 24}, "c_1", 0.5,
                  &st);
  }
}

int main(int argc, char **argv) {
  const double g = 9.81;
  const double m1 = 7.0;
  const double m2 = 5.3;
  const double m3 = 1.3;
  const double m6 = 1.5;
  const double m7 = 0.8;
  const double m8 = 0.6;
  const double m9 = 1.0;
  const double l = 0.8;
  const double a = 1.2;
  const double R = 0.24;
  const double l7 = 0.32;
  const double l8 = 0.30;
  const double kphi = 40.0;
  const double k1 = 700.0;
  const double k2 = 1600.0;
  const double mu2 = 1.16;
  const double hb = 0.05; // bar thickness
  const double wall = 0.03;
  const double phi0 = -0.3;
  const double th7 = -0.9;
  const double th8 = -1.5;
  const double a3 = -0.3; // slot direction in disc frame

  static lk_ref ref;
  static lk_sig sig;
  static lk_draw dr;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, -g);

  pil_vec2 P1 = {0, 0};
  pil_vec2 Lp = {l * sin(phi0), -l * cos(phi0)}; // left pin
  pil_add_ground(w, (pil_vec2){a + 0.45, 0.0}, (pil_vec2){-0.45, 0.0});

  // Bar with wall that carries spring and dashpot, one polygon
  pil_vec2 pts[] = {{Lp.x, Lp.y - hb / 2},     {Lp.x + a, Lp.y - hb / 2},
                    {Lp.x + a, Lp.y + hb / 2}, {Lp.x + wall, Lp.y + hb / 2},
                    {Lp.x + wall, Lp.y + 0.3}, {Lp.x, Lp.y + 0.3}};
  pil_vec2 cb = pil_poly_centroid(pts, 6);
  pil_body *bar = pil_add_body(w, &(pil_body_desc){.shape = pil_poly(pts, 6),
                                                   .mass = m1,
                                                   .pos = cb,
                                                   .flags = PIL_HIDE_COM,
                                                   .label = "m_1",
                                                   .label_dx = 200,
                                                   .label_dy = -22});
  pil_add_distance(w, pil_at(P1.x, P1.y), pil_on_at(bar, Lp.x, Lp.y),
                   PIL_LENGTH_INITIAL, 0);
  pil_add_distance(w, pil_at(a, 0), pil_on_at(bar, Lp.x + a, Lp.y),
                   PIL_LENGTH_INITIAL, 0);

  ref.bar = bar;
  ref.pin = (pil_vec2){Lp.x - cb.x, Lp.y - cb.y};
  ref.pivot = P1;
  ref.l = l;
  ref.kphi = kphi;
  pil_add_custom_force(w, torsion_force, &ref);

  // Disc 2 rolls on bar top
  double ytop = Lp.y + hb / 2;
  pil_vec2 dc = {Lp.x + 0.44, ytop + R};
  pil_body *disc = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(R),
                                                    .mass = m2,
                                                    .pos = dc,
                                                    .flags = PIL_SPOKE,
                                                    .label = "m_2",
                                                    .label_dx = -34,
                                                    .label_dy = 40});
  pil_add_rolling_on(w, disc, (pil_vec2){0, 0}, R, pil_on_at(bar, Lp.x, ytop),
                     (pil_vec2){1, 0});
  pil_element *sp1 = pil_add_spring(w, pil_on_at(bar, Lp.x + wall, dc.y),
                                    pil_on(disc, 0, 0), k1, 0.30);
  pil_element_set_label_at(sp1, "k_1", 0, -18);

  // Block 3 in a slot of disc, spring k_2 to centre
  pil_vec2 u3 = {cos(a3), sin(a3)};
  pil_body *blk = pil_add_body(
      w, &(pil_body_desc){.shape = pil_rect(0.1, 0.08),
                          .mass = m3,
                          .pos = {dc.x + 0.13 * u3.x, dc.y + 0.13 * u3.y},
                          .angle = a3,
                          .flags = PIL_HIDE_COM,
                          .label = "m_3",
                          .label_dx = 22,
                          .label_dy = -28});
  pil_add_prismatic(w, pil_on(blk, 0, 0), pil_on(disc, 0, 0), u3);
  pil_element *sp2 =
      pil_add_spring(w, pil_on(disc, 0, 0), pil_on(blk, -0.05, 0), k2, 0.06);
  pil_element_set_label_at(sp2, "k_2", 0, -16);

  // Block 6 on bar, friction mu_2
  pil_body *b6 =
      pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.12, 0.1),
                                       .mass = m6,
                                       .pos = {Lp.x + 1.0, ytop + 0.05},
                                       .flags = PIL_HIDE_COM,
                                       .label = "m_6",
                                       .label_dx = 0,
                                       .label_dy = 40});
  pil_add_floor_contact(w, b6, pil_on_at(bar, Lp.x, ytop), (pil_vec2){1, 0},
                        1e5, 400.0, mu2);

  // Double pendulum under middle of bar
  pil_vec2 S1 = {Lp.x + a / 2, Lp.y - hb / 2};
  pil_vec2 d7 = {sin(th7), -cos(th7)};
  pil_vec2 d8 = {sin(th8), -cos(th8)};
  pil_vec2 H = {S1.x + l7 * d7.x, S1.y + l7 * d7.y};
  pil_vec2 B = {H.x + l8 * d8.x, H.y + l8 * d8.y};
  pil_body *r7 = pil_add_body(
      w,
      &(pil_body_desc){.shape = pil_rect(l7, 0.03),
                       .mass = m7,
                       .pos = {S1.x + 0.5 * l7 * d7.x, S1.y + 0.5 * l7 * d7.y},
                       .angle = th7 - M_PI / 2,
                       .flags = PIL_HIDE_COM,
                       .fill = 0.3});
  pil_body *r8 = pil_add_body(
      w, &(pil_body_desc){.shape = pil_rect(l8, 0.03),
                          .mass = m8,
                          .pos = {H.x + 0.5 * l8 * d8.x, H.y + 0.5 * l8 * d8.y},
                          .angle = th8 - M_PI / 2,
                          .flags = PIL_HIDE_COM,
                          .fill = 0.3});
  pil_body *ball = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(0.045),
                                                    .mass = m9,
                                                    .pos = B,
                                                    .flags = PIL_HIDE_COM,
                                                    .label = "m_9",
                                                    .label_dx = -4,
                                                    .label_dy = 34});
  pil_add_revolute(w, pil_on_at(bar, S1.x, S1.y), pil_on_at(r7, S1.x, S1.y));
  pil_add_revolute(w, pil_on_at(r7, H.x, H.y), pil_on_at(r8, H.x, H.y));
  pil_add_weld(w, pil_on_at(r8, B.x, B.y), pil_on_at(ball, B.x, B.y));

  // Annotations
  pil_add_line(w, pil_at(P1.x, P1.y), pil_at(P1.x, -l - 0.1), 1);
  pil_add_angle_between(w, pil_at(P1.x, P1.y), pil_at(P1.x, -0.5),
                        pil_on_at(bar, Lp.x, Lp.y), 60, "\xcf\x86");
  pil_add_dimension(w, pil_at(P1.x, P1.y), pil_on_at(bar, Lp.x, Lp.y), 24, "l");
  pil_add_angle_marker(w, pil_on_at(bar, S1.x, S1.y), -M_PI / 2, r7, 34,
                       "\xce\xb8_7");
  pil_add_angle_marker(w, pil_on_at(r7, H.x, H.y), -M_PI / 2, r8, 30,
                       "\xce\xb8_8");
  pil_add_dimension(w, pil_on_at(bar, Lp.x, Lp.y - 0.5),
                    pil_on_at(bar, S1.x, Lp.y - 0.5), 0, "a/2");
  pil_add_dimension(w, pil_on_at(bar, S1.x, Lp.y - 0.5),
                    pil_on_at(bar, Lp.x + a, Lp.y - 0.5), 0, "a/2");
  pil_add_gravity_indicator(w, (pil_vec2){-0.35, -0.35});

  sig.bar = bar;
  sig.disc = disc;
  sig.blk = blk;
  sig.r7 = r7;
  sig.r8 = r8;
  sig.pin = ref.pin;
  sig.pivot = P1;
  sig.slot = u3;
  sig.wall = wall;
  dr.s = &sig;
  dr.a = a;
  dr.l = l;
  pil_world_add_overlay(w, lk_overlay, &dr);

  pil_plot *pa =
      pil_add_plot(w, 0.68, 0.05, 0.99, 0.47, "t [s]", "angle [rad]", 6.0);
  pil_plot_add_signal(pa, "\xcf\x86", sig_phi, &sig);
  pil_plot_add_signal(pa, "\xce\xb8_7", sig_th7, &sig);
  pil_plot_add_signal(pa, "\xce\xb8_8", sig_th8, &sig);
  pil_plot *pb =
      pil_add_plot(w, 0.68, 0.55, 0.99, 0.97, "t [s]", "position [m]", 6.0);
  pil_plot_add_signal(pb, "disc", sig_disc, &sig);
  pil_plot_add_signal(pb, "slider", sig_slot, &sig);
  pil_plots_record(w, 0.03, 5e-4);

  pil_view view;
  pil_view_fit_region(&view, -0.6, -1.9, a + 0.5, 0.15, 1500, 800, 0.0, 0.0,
                      0.66, 1.0, 10);
  int rc = example_main(w, &view, "Parallelogram linkage with a rolling disc",
                        argc, argv, "linkage.png");

  pil_world_free(w);

  return rc;
}
