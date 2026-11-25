/*
 * shaker.c - plan view: box on springs, a rolling disc and a pendulum
 * Gravity is normal to page, so nothing is pulled down: all motion comes from
 * springs and from prescribed displacement u(t) of right support
 * Box 3m is held by spring k and dashpot d to wall and by spring 2k to driven
 * point; on it a disc 2m rolls against spring k and dashpot d; a pendulum m of
 * length l hangs from ceiling hinge and is tied to box by spring 2k
 * Three coordinates: box x, disc position, pendulum angle phi
 * Everything is public-API
 */

#include "example_main.h"
#include <math.h>

typedef struct {
  pil_body *box;
  pil_body *disc;
  pil_body *bob;
  pil_vec2 hinge;
  double u0, U, Om;
} sh_sig;

static double drive_u(double t, void *user) {
  const sh_sig *s = user;

  return s->u0 + s->U * sin(s->Om * t);
}

static double sig_u(const pil_world *w, void *user) {
  const sh_sig *s = user;

  return drive_u(pil_world_time(w), user) - s->u0;
}

static double sig_box(const pil_world *w, void *user) {
  (void)w;

  return ((const sh_sig *)user)->box->pos.x -
         ((const sh_sig *)user)->box->pos0.x;
}

static double sig_disc(const pil_world *w, void *user) {
  (void)w;

  return ((const sh_sig *)user)->disc->pos.x -
         ((const sh_sig *)user)->disc->pos0.x;
}

static double sig_phi(const pil_world *w, void *user) {
  const sh_sig *s = user;
  (void)w;

  return atan2(s->bob->pos.x - s->hinge.x, s->hinge.y - s->bob->pos.y);
}

int main(int argc, char **argv) {
  const double m = 1.0;
  const double k = 80.0;
  const double d = 5.0;
  const double R = 0.16;
  const double l = 0.7;
  const double phi0 = 0.3;
  static sh_sig sig;
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);

  pil_add_ground(w, (pil_vec2){1.7, 1.5}, (pil_vec2){-1.3, 1.5});
  pil_add_ground(w, (pil_vec2){-1.3, 1.5}, (pil_vec2){-1.3, -0.35});
  pil_add_gravity_into_page(w, (pil_vec2){1.35, 1.35});

  // Box with two posts on top, one polygon
  pil_vec2 pts[] = {{-0.5, 0},     {0.5, 0},      {0.5, 0.5},   {0.47, 0.5},
                    {0.47, 0.66},  {0.43, 0.66},  {0.43, 0.5},  {-0.43, 0.5},
                    {-0.43, 0.66}, {-0.47, 0.66}, {-0.47, 0.5}, {-0.5, 0.5}};
  pil_vec2 cb = pil_poly_centroid(pts, 12);
  pil_body *box = pil_add_body(w, &(pil_body_desc){.shape = pil_poly(pts, 12),
                                                   .mass = 3 * m,
                                                   .pos = cb,
                                                   .flags = PIL_HIDE_COM,
                                                   .label = "3m",
                                                   .label_dy = -10});
  pil_add_prismatic(w, pil_on(box, 0, 0), pil_at(cb.x, cb.y), (pil_vec2){1, 0});

  pil_element *s1 = pil_add_spring(
      w, pil_at(-1.3, 0.35), pil_on_at(box, -0.5, 0.35), k, PIL_REST_INITIAL);
  pil_element_set_label_at(s1, "k", 0, -16);
  pil_element *d1 =
      pil_add_damper(w, pil_at(-1.3, 0.1), pil_on_at(box, -0.5, 0.1), d);
  pil_element_set_label_at(d1, "d", 0, -16);

  // Driven support: a point forced to u(t) along a horizontal guide
  sig.u0 = 1.3;
  sig.U = 0.1;
  sig.Om = 5.0;
  pil_body *drv =
      pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.025),
                                       .mass = 0.1,
                                       .pos = {sig.u0, 0.25},
                                       .vel = {sig.U * sig.Om, 0},
                                       .flags = PIL_HIDE_COM});
  pil_add_prismatic(w, pil_on(drv, 0, 0), pil_at(sig.u0, 0.25),
                    (pil_vec2){1, 0});
  pil_add_coord_drive(w, drv, PIL_COORD_X, drive_u, &sig);
  pil_element *s2 = pil_add_spring(w, pil_on_at(box, 0.5, 0.25),
                                   pil_on(drv, 0, 0), 2 * k, PIL_REST_INITIAL);
  pil_element_set_label_at(s2, "2k", 0, -16);
  pil_add_arrow(w, pil_on(drv, 0, 0), 40, 0, "u(t)");

  // Disc on box, rolling without slipping
  pil_body *disc = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(R),
                                                    .mass = 2 * m,
                                                    .pos = {-0.15, 0.5 + R},
                                                    .flags = PIL_SPOKE,
                                                    .label = "2m",
                                                    .label_dx = 0,
                                                    .label_dy = 28});
  pil_add_rolling_on(w, disc, (pil_vec2){0, 0}, R, pil_on_at(box, -0.5, 0.5),
                     (pil_vec2){1, 0});
  pil_element *s3 = pil_add_spring(w, pil_on_at(box, -0.43, 0.5 + R),
                                   pil_on(disc, 0, 0), k, PIL_REST_INITIAL);
  pil_element_set_label_at(s3, "k", 0, 14);
  pil_element *d3 =
      pil_add_damper(w, pil_on(disc, 0, 0), pil_on_at(box, 0.43, 0.5 + R), d);
  pil_element_set_label_at(d3, "d", 0, 14);

  // Pendulum from ceiling, tied to box by spring 2k
  pil_vec2 H = {0.15, 1.5};
  pil_vec2 B = {H.x + l * sin(phi0), H.y - l * cos(phi0)};
  pil_body *bob =
      pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.05),
                                       .mass = m,
                                       .pos = B,
                                       .label = "m",
                                       .label_dx = 16,
                                       .label_dy = 10});
  pil_add_distance(w, pil_at(H.x, H.y), pil_on(bob, 0, 0), PIL_LENGTH_INITIAL,
                   0);
  pil_element *s4 = pil_add_spring(w, pil_on(bob, 0, 0),
                                   pil_on_at(box, 0.43, 0.5 + R), 2 * k, 0.25);
  pil_element_set_label_at(s4, "2k", 24, 0);
  pil_add_line(w, pil_at(H.x, H.y), pil_at(H.x, H.y - l - 0.1), 1);
  pil_add_angle_between(w, pil_at(H.x, H.y), pil_at(H.x, H.y - 0.4),
                        pil_on(bob, 0, 0), 60, "\xcf\x86");
  pil_add_dimension(w, pil_at(H.x, H.y), pil_on(bob, 0, 0), 22, "l");

  sig.box = box;
  sig.disc = disc;
  sig.bob = bob;
  sig.hinge = H;

  pil_plot *p1 =
      pil_add_plot(w, 0.66, 0.05, 0.99, 0.47, "t [s]", "displacement [m]", 8.0);
  pil_plot_add_signal(p1, "u", sig_u, &sig);
  pil_plot_add_signal(p1, "box", sig_box, &sig);
  pil_plot_add_signal(p1, "disc", sig_disc, &sig);
  pil_plot *p2 =
      pil_add_plot(w, 0.66, 0.55, 0.99, 0.97, "t [s]", "angle [rad]", 8.0);
  pil_plot_add_signal(p2, "\xcf\x86", sig_phi, &sig);
  pil_plots_record(w, 0.03, 5e-4);

  pil_view view;
  pil_view_fit_region(&view, -1.45, -0.45, 1.65, 1.65, 1500, 700, 0.0, 0.0,
                      0.64, 1.0, 10);
  int rc = example_main(w, &view, "Driven box with rolling disc and pendulum",
                        argc, argv, "shaker.png");

  pil_world_free(w);

  return rc;
}
