/*
 * atwood.c - two masses over a pulley, with force arrows and a floor stop
 */
#include "example_main.h"
#include <math.h>

static void pulley_pin(cairo_t *cr, const pil_view *v, void *user) {
  pil_style st = pil_style_textbook();
  const pil_vec2 *c = user;
  pil_draw_hinge(cr, pil_view_to_screen(v, *c), &st);
}

int main(int argc, char **argv) {
  const double R = 0.2, yc = 2.0, h = 0.28;
  static pil_vec2 centre = {0.0, 2.0};
  pil_world *w = pil_world_new();

  /* ceiling, hanger and pulley */
  pil_add_ground(w, (pil_vec2){0.9, yc + 0.7}, (pil_vec2){-0.9, yc + 0.7});
  pil_add_body(w, &(pil_body_desc){.shape = pil_rect(0.05, 0.7),
                                   .mass = 1,
                                   .pos = {0, yc + 0.35},
                                   .flags = pil_STATIC | pil_HIDE_COM});
  pil_add_body(w, &(pil_body_desc){.shape = pil_disc(R),
                                   .mass = 1,
                                   .pos = centre,
                                   .flags = pil_STATIC | pil_HIDE_COM});
  pil_world_add_overlay(w, pulley_pin, &centre);
  pil_add_ground(w, (pil_vec2){-0.9, 0}, (pil_vec2){0.9, 0});

  // two masses, each guided on a vertical line, tied by a rope over pulley
  pil_body *A = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(h, h),
                                                 .mass = 3.0,
                                                 .pos = {-R, 0.95},
                                                 .label = "m_1",
                                                 .label_dx = -40,
                                                 .label_dy = 0});
  pil_body *B = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(h, h),
                                                 .mass = 2.0,
                                                 .pos = {+R, 0.62},
                                                 .label = "m_2",
                                                 .label_dx = 40,
                                                 .label_dy = 0});
  pil_add_prismatic(w, pil_on(A, 0, 0), pil_at(-R, 0), (pil_vec2){0, 1});
  pil_add_prismatic(w, pil_on(B, 0, 0), pil_at(+R, 0), (pil_vec2){0, 1});
  pil_add_floor_contact(w, A, pil_at(0, 0), (pil_vec2){1, 0}, 2e5, 2e3, 0.0);
  pil_add_floor_contact(w, B, pil_at(0, 0), (pil_vec2){1, 0}, 2e5, 2e3, 0.0);
  pil_add_pulley_rope(w, pil_on(A, 0, h / 2), pil_on(B, 0, h / 2),
                      (pil_vec2){-R, yc}, (pil_vec2){+R, yc},
                      pil_LENGTH_INITIAL);

  // free-body style arrows
  pil_add_arrow(w, pil_on(A, 0, -h / 2), 0, -46, "m_1g");
  pil_add_arrow(w, pil_on(B, 0, -h / 2), 0, -46, "m_2g");
  pil_add_gravity_indicator(w, (pil_vec2){1.05, 2.4});

  pil_view view;
  pil_view_fit(&view, -1.1, -0.45, 1.3, 2.8, 520, 640, 22);
  int rc = example_main(w, &view, "Atwood machine", argc, argv, "atwood.png");
  pil_world_free(w);

  return rc;
}
