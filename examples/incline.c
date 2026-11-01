/*
 * incline.c - a "roof" wedge
 * Left slope: a cylinder rolls without slipping while tied to apex post bya
 * spring
 * Right slope: a block with Coulomb friction on same kind of spring
 * Everything is public-API
 */

#include "example_main.h"
#include <math.h>

int main(int argc, char **argv) {
  const double base = 3.2;
  const double rise = 1.0;
  const pil_vec2 apex = {base / 2, rise};
  const double al = atan2(rise, base / 2);
  const pil_vec2 ul = {cos(al), sin(al)};  // left slope: up-hill
  const pil_vec2 nl = {-sin(al), cos(al)}; // left slope: normal
  const pil_vec2 ur = {cos(al), -sin(al)}; // right slope: down-hill
  const pil_vec2 nr = {sin(al), cos(al)};  // right slope: normal

  pil_world *w = pil_world_new();

  // wedge (static polygon; body origin = centroid) and apex post
  pil_vec2 tri[3] = {{0, 0}, {base, 0}, apex};
  pil_add_body(w, &(pil_body_desc){.shape = pil_poly(tri, 3),
                                   .mass = 1,
                                   .pos = {base / 2, rise / 3},
                                   .flags = PIL_STATIC | PIL_HIDE_COM});
  pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.24, 0.5),
                                   .mass = 1,
                                   .pos = {apex.x, apex.y + 0.05},
                                   .flags = PIL_STATIC | PIL_HIDE_COM});
  pil_add_ground(w, (pil_vec2){-0.5, 0}, (pil_vec2){base + 0.5, 0});

  // left: cylinder rolling on slope, spring to post
  const double r = 0.15;
  pil_vec2 cpos = {apex.x - ul.x * 1.0 + nl.x * r,
                   apex.y - ul.y * 1.0 + nl.y * r};
  pil_body *cyl = pil_add_body(w, &(pil_body_desc){.shape = pil_disc(r),
                                                   .mass = 3.0,
                                                   .pos = cpos,
                                                   .flags = PIL_SPOKE,
                                                   .label = "M",
                                                   .label_dx = -18,
                                                   .label_dy = 18});
  pil_add_rolling(w, cyl, (pil_vec2){0, 0}, r, (pil_vec2){0, 0}, ul);
  pil_vec2 ca = {apex.x + nl.x * r, apex.y + nl.y * r};
  pil_element *s1 =
      pil_add_spring(w, pil_at(ca.x, ca.y), pil_on(cyl, 0, 0), 60.0, 0.55);
  pil_element_set_label_at(s1, "k_1", -14, 24);

  // right: block with friction on slope, spring to post
  const double hb = 0.13;
  pil_vec2 bpos = {apex.x + ur.x * 1.0 + nr.x * hb,
                   apex.y + ur.y * 1.0 + nr.y * hb};
  pil_body *blk =
      pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.42, 2 * hb),
                                       .mass = 2.0,
                                       .pos = bpos,
                                       .angle = -al,
                                       .label = "m",
                                       .label_dx = 22,
                                       .label_dy = 22});
  pil_add_floor_contact(w, blk, pil_at(apex.x, apex.y), (pil_vec2){ur.x, ur.y},
                        3e5, 2e3, 0.2);
  pil_vec2 ba = {apex.x + nr.x * hb, apex.y + nr.y * hb};
  pil_element *s2 =
      pil_add_spring(w, pil_at(ba.x, ba.y), pil_on(blk, -0.21, 0), 40.0, 0.5);
  pil_element_set_label_at(s2, "k_2", 14, 26);

  pil_add_label(w, pil_at(0.62, 0.10), "\xce\xb1", 0, 0);
  pil_add_label(w, pil_at(base - 0.62, 0.10), "\xce\xb1", 0, 0);
  pil_add_label(w, pil_on(blk, 0, 0), "\xce\xbc", 0, -34);
  pil_add_gravity_indicator(w, (pil_vec2){base + 0.25, 1.55});

  pil_view view;
  pil_view_fit(&view, -0.55, -0.15, base + 0.55, 1.75, 900, 470, 24);
  int rc = example_main(w, &view, "Cylinder and block on a roof wedge", argc,
                        argv, "incline.png");
  pil_world_free(w);

  return rc;
}
