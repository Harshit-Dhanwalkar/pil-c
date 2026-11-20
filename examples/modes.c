/*
 * modes.c - two coupled oscillators and their normal modes
 * Three identical chains wall - k - m - k_c - m - k - wall
 * Top: both masses start displaced together and stay in step (omega_1 = sqrt(k
 * / m)); middle: displaced in opposition (omega_2 = sqrt((k + 2 k_c) / m));
 * bottom: only left mass is displaced, a sum of both modes, so energy beats
 * between masses
 * Plot: two displacements of bottom chain
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 3.0
#include "example_main.h"

typedef struct {
  pil_body *b;
  double eq;
} probe;

static double sig_probe(const pil_world *w, void *user) {
  (void)w;
  const probe *p = user;

  return p->b->pos.x - p->eq;
}

static const double m = 1.0, k = 4.0, kc = 2.0;
static const double L = 1.2, hw = 0.2, A = 0.25;

static void chain(pil_world *w, double y, double a1, double a2, probe *pr) {
  pil_body *b[2];
  for (int i = 0; i < 2; i++) {
    double eq = (i + 1) * L;
    b[i] = pil_add_body(w, &(pil_body_desc){.shape = pil_rect(2 * hw, 0.3),
                                            .mass = m,
                                            .pos = {eq + (i ? a2 : a1), y},
                                            .fill = i ? 0.6 : 0.85,
                                            .flags = PIL_HIDE_COM});
    pil_add_prismatic(w, pil_on(b[i], 0, 0), pil_at(eq, y), (pil_vec2){1, 0});
    if (pr) {
      pr[i].b = b[i];
      pr[i].eq = eq;
    }
  }

  pil_add_ground(w, (pil_vec2){0, y + 0.15}, (pil_vec2){0, y - 0.35});
  pil_add_ground(w, (pil_vec2){3 * L, y + 0.15}, (pil_vec2){3 * L, y - 0.35});
  pil_add_ground(w, (pil_vec2){-0.1, y - 0.15},
                 (pil_vec2){3 * L + 0.1, y - 0.15});
  pil_add_spring(w, pil_at(0, y), pil_on(b[0], -hw, 0), k, L - hw);
  pil_add_spring(w, pil_on(b[0], hw, 0), pil_on(b[1], -hw, 0), kc, L - 2 * hw);
  pil_add_spring(w, pil_on(b[1], hw, 0), pil_at(3 * L, y), k, L - hw);
}

int main(int argc, char **argv) {
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  static probe pr[2];

  chain(w, 3.0, A, A, NULL);
  chain(w, 1.6, A, -A, NULL);
  chain(w, 0.2, A, 0, pr);
  pil_add_label(w, pil_at(0.0, 3.55),
                "mode 1: in phase, \xcf\x89 = sqrt(k / m)", 2, 0);
  pil_add_label(w, pil_at(0.0, 2.15),
                "mode 2: opposite, \xcf\x89 = sqrt((k + 2 k_c) / m)", 2, 0);
  pil_add_label(w, pil_at(0.0, 0.75), "sum of both modes: energy beats", 2, 0);

  pil_plot *pl = pil_add_plot(w, 0.60, 0.15, 0.98, 0.85, "t [s]",
                              "displacement [m]", 12.0);
  pil_plot_add_signal(pl, "x_1", sig_probe, &pr[0]);
  pil_plot_add_signal(pl, "x_2", sig_probe, &pr[1]);
  pil_plots_record(w, 0.05, 5e-4);

  pil_view view;
  pil_view_fit_region(&view, -0.3, -0.5, 3.9, 3.9, 1500, 700, 0.0, 0.0, 0.58,
                      1.0, 12);
  int rc =
      example_main(w, &view, "Coupled oscillators", argc, argv, "modes.png");

  pil_world_free(w);

  return rc;
}
