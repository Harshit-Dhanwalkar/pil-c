/*
 * halfcar.c - "half-car on a wavy road" figure, built only from public API of
 * library
 *
 *   ./halfcar [frame.png] | --record out.gif [seconds] | --live
 */

#include "example_main.h"
#include <math.h>

#define PX_PER_M 215.0

// road profile h(s), s = x + v0 * t
static double road(double s, void *user) {
  (void)user;
  return 0.045 * sin(2 * M_PI * s / 2.1) +
         0.012 * sin(2 * M_PI * s / 0.65 + 1.0);
}

// Anything library does not provide can be drawn in an overlay:
// here dash-dot axis through COM, its pin, and v_0 arrow
static void body_axis_overlay(cairo_t *cr, const pil_view *v, void *user) {
  const pil_body *M = user;
  pil_style st = pil_style_textbook();
  pil_vec2 l = pil_view_to_screen(v, pil_v2(-1.78, M->pos.y));
  pil_vec2 r = pil_view_to_screen(v, pil_v2(1.35, M->pos.y));
  cairo_save(cr);
  cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);
  cairo_set_line_width(cr, 1.0);
  double dash[4] = {14, 3.5, 2, 3.5};
  cairo_set_dash(cr, dash, 4, 0);
  cairo_move_to(cr, l.x, l.y);
  cairo_line_to(cr, r.x, r.y);
  cairo_stroke(cr);
  cairo_set_dash(cr, NULL, 0, 0);
  pil_draw_hinge(cr, l, &st);
  pil_draw_arrow(cr, r, pil_v2(r.x + 52, r.y), 10, &st);
  pil_draw_text(cr, pil_v2(r.x + 62, r.y + 2), "v_0", 0.0, &st);
  cairo_restore(cr);
}

int main(int argc, char **argv) {
  // geometry [m]
  const double b = 2.4;
  const double a = 0.36;
  const double l = 0.9;
  const double r = 0.2;
  // speed [m/s]
  const double v0 = 5.0;
  // tyre
  const double kt = 2.2e5;
  const double ct = 400;

  pil_world *w = pil_world_new();
  pil_terrain *ground =
      pil_add_terrain(w, road, NULL, v0); // road scrolls under car

  // sprung mass M: free to heave and pitch (slot on a vertical line through
  // COM)
  pil_body *M =
      pil_add_body(w, &(pil_body_desc){.shape = pil_rect(b, a),
                                       .mass = 1000,
                                       .inertia = 1100,
                                       .pos = {0, r + 0.02 + 0.62 + a / 2},
                                       .label = "M",
                                       .label_dx = 40,
                                       .label_dy = 0});
  pil_add_slot(w, pil_on(M, 0, 0), pil_at(0, 0), (pil_vec2){0, 1});

  // two axles; each = hidden non-rotating carrier (spring/damper end) +
  // spinning wheel
  const double xs[2] = {-l, l};
  const double mw[2] = {40, 45};
  const double ks[2] = {32000, 38000};
  const double ds[2] = {2400, 2800};
  const char *wl[2] = {"m_1", "m_2"};
  const double *kl[2] = {"k_1", "k_2"};
  const double *dl[2] = {"d_1", "d_2"};
  for (int i = 0; i < 2; i++) {
    pil_vec2 p = {xs[i], r + road(xs[i], NULL)};
    pil_body *car =
        pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.01),
                                         .mass = 0.5,
                                         .pos = p,
                                         .flags = PIL_HIDDEN});
    pil_body *wh =
        pil_add_body(w, &(pil_body_desc){.shape = pil_disc(r),
                                         .mass = mw[i],
                                         .pos = p,
                                         .flags = PIL_SPOKE | PIL_HIDE_COM,
                                         .label = wl[i],
                                         .label_dx = 0,
                                         .label_dy = -(r * PX_PER_M + 20)});
    pil_add_prismatic(w, pil_on(car, 0, 0), pil_at(xs[i], 0), (pil_vec2){0, 1});
    pil_element_hide(pil_add_revolute(w, pil_on(car, 0, 0), pil_on(wh, 0, 0)));
    pil_add_rotation_drive(w, wh, -v0 / r); // rolling forward
    pil_add_terrain_contact(w, ground, wh, r, kt, ct);

    pil_element *sp =
        pil_add_spring(w, pil_on(M, xs[i] - 0.09, -a / 2),
                       pil_on(car, -0.09, r + 0.02), ks[i], PIL_REST_INITIAL);
    pil_element *dp = pil_add_damper(w, pil_on(M, xs[i] + 0.09, -a / 2),
                                     pil_on(car, 0.09, r + 0.02), ds[i]);
    pil_element_set_label(sp, kl[i]);
    pil_element_set_label_at(dp, dl[i], 32, 0);
  }

  // textbook annotations
  pil_add_dimension(w, pil_on(M, -b / 2, a / 2), pil_on(M, b / 2, a / 2), 38,
                    "b");
  pil_add_dimension(w, pil_on(M, -l, a / 2), pil_on(M, 0, a / 2), 16,
                    "\xe2\x84\x93");
  pil_add_dimension(w, pil_on(M, 0, a / 2), pil_on(M, l, a / 2), 16,
                    "\xe2\x84\x93");
  pil_add_dimension(w, pil_on(M, -b / 2, -a / 2), pil_on(M, -b / 2, a / 2), 24,
                    "a");
  pil_add_angle_marker(w, pil_on(M, 0, 0), 0.0, M, 190, "\xcf\x86");
  pil_add_gravity_indicator(w, pil_v2(1.8, 1.45));
  pil_world_add_overlay(w, body_axis_overlay, M);

  pil_world_settle(w, 6.0); /* start from static equilibrium */

  pil_view view;
  pil_view_fit(&view, -1.95, -0.12, 2.05, 1.62, 900, 440, 24);
  int rc = example_main(w, &view, "Half-car on a wavy road", argc, argv,
                        "halfcar.png");
  pil_world_free(w);
  return rc;
}
