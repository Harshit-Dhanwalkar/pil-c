/*
 * bowl.c - a bowl (annular sector) rolling on ground between two springs, a
 * ball rolling inside it, and a T-shaped pendulum with two end masses hinged at
 * ball's centre and tied to bowl by a torsion spring
 * Public API only
 */

#include "example_main.h"
#include <math.h>

int main(int argc, char **argv) {
  const double Ro = 0.62, Ri = 0.48; // bowl: outer/inner radius
  const double beta = 0.30;          // end angle above horizontal
  const double r = 0.09;             // ball radius
  const double D = Ri - r;           // centre distance from bowl centre
  const double ds = 1.45, wp = 0.26; // spring anchors, bar width,
  const double Lrod = 0.58;          //  rod length
  const double Rm = 0.5 * (Ri + Ro);
  const double x0 = 0.10;               // initial bowl offset
  const double psi0 = -M_PI / 2 + 0.55; // initial ball angle */

  pil_world *w = pil_world_new();
  pil_add_ground(w, (pil_vec2){-1.85, 0}, (pil_vec2){1.85, 0});

  // bowl: arc band from (-\pi - \beta) to \beta, centroid computed
  pil_body *bowl = pil_add_body(
      w, &(pil_body_desc){.shape = pil_arc_band(Ri, Ro, -M_PI - beta, beta),
                          .mass = 4.0});
  pil_vec2 o = pil_body_shape_origin(bowl); // circle centre in body frame
  pil_body_set_state(bowl, (pil_vec2){x0 - o.x, Ro - o.y}, 0, (pil_vec2){0, 0},
                     0);
  pil_add_rolling(w, bowl, o, Ro, (pil_vec2){0, 0},
                  (pil_vec2){1, 0}); // rolls on ground

  const double k = 90.0; // springs from bowl's ends to ground
  double rest = hypot(ds - Rm * cos(beta), Ro + Rm * sin(beta));
  for (int s = -1; s <= 1; s += 2) {
    pil_element *sp = pil_add_spring(
        w, pil_at(s * ds, 0),
        pil_on(bowl, o.x + s * Rm * cos(beta), o.y + Rm * sin(beta)), k, rest);

    pil_element_set_label_at(sp, "k", s * 28, 16);
  }

  // ball, rolling on inner surface
  pil_vec2 C0 = {x0 + D * cos(psi0), Ro + D * sin(psi0)};
  pil_body *ball =
      pil_add_body(w, &(pil_body_desc){.shape = pil_disc(r),
                                       .mass = 1.5,
                                       .pos = C0,
                                       .flags = PIL_SPOKE | PIL_HIDE_COM,
                                       .label = "m_r",
                                       .label_dx = 16,
                                       .label_dy = 18});

  pil_add_circle_rolling(w, pil_on(ball, 0, 0), r, pil_on(bowl, o.x, o.y), Ri,
                         1);

  // pendulum: a T-shaped rod (polygon, hinge at shape origin) + two welded
  // masses
  const double hw = 0.015;
  const double th = 0.03;
  pil_vec2 T[8] = {{-hw, 0},
                   {-hw, -Lrod},
                   {-wp / 2, -Lrod},
                   {-wp / 2, -Lrod - th},
                   {wp / 2, -Lrod - th},
                   {wp / 2, -Lrod},
                   {hw, -Lrod},
                   {hw, 0}};
  pil_body *rod = pil_add_body(w, &(pil_body_desc){.shape = pil_poly(T, 8),
                                                   .mass = 0.3,
                                                   .flags = PIL_HIDE_COM});
  pil_vec2 ro = pil_body_shape_origin(rod); // hinge, in rod's body frame
  pil_body_set_state(rod, (pil_vec2){C0.x - ro.x, C0.y - ro.y}, 0,
                     (pil_vec2){0, 0}, 0);
  pil_add_revolute(w, pil_on(ball, 0, 0), pil_on(rod, ro.x, ro.y));
  for (int s = -1; s <= 1; s += 2) {
    pil_vec2 off = {s * wp / 2, -Lrod - 0.5 * th};
    pil_body *mp =
        pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.04),
                                         .mass = 0.8,
                                         .pos = {C0.x + off.x, C0.y + off.y},
                                         .label = "m_p",
                                         .label_dx = s * 20,
                                         .label_dy = -14});
    pil_add_weld(w, pil_on(rod, ro.x + off.x, ro.y + off.y), pil_on(mp, 0, 0));
  }

  pil_add_torsion_spring(w, rod, bowl, 1.2, 0.0, 0.08);
  pil_add_label(w, pil_on(ball, 0, 0), "k_\xcf\x86", -24, -10);

  // annotations
  pil_add_dimension(w, pil_at(0, 0), pil_at(-ds, 0), 150, "d_s");
  pil_add_dimension(w, pil_at(0, 0), pil_at(ds, 0), -150, "d_s");
  pil_add_dimension(w, pil_on(rod, ro.x - wp / 2, ro.y - Lrod - th - 0.02),
                    pil_on(rod, ro.x + wp / 2, ro.y - Lrod - th - 0.02), -26,
                    "w_p");
  const double a_ro = M_PI + 0.28;
  const double a_ri = M_PI + 0.80; // radius dimensions end on bowl's surfaces
  pil_add_dimension(w, pil_on(bowl, o.x, o.y),
                    pil_on(bowl, o.x + Ro * cos(a_ro), o.y + Ro * sin(a_ro)), 0,
                    "R_o");
  pil_add_dimension(w, pil_on(bowl, o.x, o.y),
                    pil_on(bowl, o.x + Ri * cos(a_ri), o.y + Ri * sin(a_ri)), 0,
                    "R_i");
  pil_add_line(w, pil_on(bowl, o.x, o.y), pil_on(bowl, o.x + 1.12 * Ro, o.y),
               1);
  pil_add_angle_between(
      w, pil_on(bowl, o.x, o.y), pil_on(bowl, o.x + Ro, o.y),
      pil_on(bowl, o.x + Ro * cos(beta), o.y + Ro * sin(beta)), 105,
      "\xce\xb2_s");
  pil_add_leader(w, pil_on(bowl, o.x + Rm * cos(-0.7), o.y + Rm * sin(-0.7)),
                 62, -4, "M");
  pil_add_gravity_indicator(w, (pil_vec2){1.55, 1.22});

  pil_view view;
  pil_view_fit(&view, -1.8, -0.80, 1.8, 1.38, 900, 560, 24);
  int rc = example_main(w, &view, "Rolling bowl with ball and pendulum", argc,
                        argv, "bowl.png");
  pil_world_free(w);
  return rc;
}
