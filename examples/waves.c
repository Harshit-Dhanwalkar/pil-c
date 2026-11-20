/*
 * waves.c - standing waves on a chain of beads (normal modes of many masses)
 * N beads of mass m, joined by springs of zero rest length (so transverse
 * restoring force is exactly linear in displacements), ends fixed
 * Mode n has shape sin(n pi j / (N + 1)) and frequency omega_n = 2 sqrt(k / m)
 * sin(n pi / (2 (N + 1)))
 * Three chains start in modes 1, 2 and 3: each bead moves in step, nodes stay
 * still and higher mode faster oscillation
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 0.4
#define PIL_EXAMPLE_IMAGE_TIME 4.0
#include "example_main.h"
#include <math.h>

#define NB 19

static const double K = 100.0;
static const double M = 1.0;
static const double LEN = 8.0;
static const double AMP = 0.4;

typedef struct {
  pil_body *b[NB];
  double y0;
  int mode;
} chain_t;

static chain_t chains[3];

static void draw_chains(cairo_t *cr, const pil_view *v, void *user) {
  (void)user;
  pil_style st = pil_style_textbook();
  cairo_save(cr);

  for (int c = 0; c < 3; c++) {
    // Amplitude envelopes
    cairo_set_source_rgb(cr, 0.7, 0.7, 0.7);
    cairo_set_line_width(cr, st.thin_w * 0.8);
    double dash[] = {4.0, 4.0};
    cairo_set_dash(cr, dash, 2, 0);
    for (int sign = -1; sign <= 1; sign += 2) {
      cairo_new_path(cr);
      int pts = 100;
      for (int i = 0; i <= pts; i++) {
        double x = LEN * i / pts;
        double env = sign * AMP * sin(chains[c].mode * M_PI * x / LEN);
        pil_vec2 p = pil_view_to_screen(v, pil_v2(x, chains[c].y0 + env));
        if (i == 0)
          cairo_move_to(cr, p.x, p.y);
        else
          cairo_line_to(cr, p.x, p.y);
      }
      cairo_stroke(cr);
    }
    cairo_set_dash(cr, NULL, 0, 0);

    // String connecting beads
    cairo_set_source_rgb(cr, st.ink[0], st.ink[1], st.ink[2]);
    cairo_set_line_width(cr, st.thin_w + 0.5);
    pil_vec2 p = pil_view_to_screen(v, pil_v2(0, chains[c].y0));
    cairo_move_to(cr, p.x, p.y);
    for (int j = 0; j < NB; j++) {
      p = pil_view_to_screen(v, chains[c].b[j]->pos);
      cairo_line_to(cr, p.x, p.y);
    }
    p = pil_view_to_screen(v, pil_v2(LEN, chains[c].y0));
    cairo_line_to(cr, p.x, p.y);
    cairo_stroke(cr);

    // Mark nodes (where displacement is always 0)
    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
    for (int j = 1; j < NB + 1; j++) {
      double x = j * LEN / (NB + 1);
      // Node condition: sin(n pi x / L) == 0
      if (fabs(sin(chains[c].mode * M_PI * j / (NB + 1))) < 0.05) {
        pil_vec2 np = pil_view_to_screen(v, pil_v2(x, chains[c].y0));
        cairo_arc(cr, np.x, np.y, 4.0, 0, 2 * M_PI);
        cairo_fill(cr);
      }
    }
  }

  cairo_restore(cr);
}

static void chain(pil_world *w, double y0, int mode) {
  chain_t *ch = &chains[mode - 1];
  ch->y0 = y0;
  ch->mode = mode;
  double dx = LEN / (NB + 1);
  pil_body **b = ch->b;
  for (int j = 0; j < NB; j++) {
    double x = (j + 1) * dx;
    double dy = AMP * sin(mode * M_PI * (j + 1) / (NB + 1));
    b[j] = pil_add_body(w, &(pil_body_desc){.shape = pil_point_mass(0.07),
                                            .mass = M,
                                            .pos = {x, y0 + dy},
                                            .flags = PIL_HIDE_COM});
    pil_add_prismatic(w, pil_on(b[j], 0, 0), pil_at(x, y0), (pil_vec2){0, 1});
  }

  // Trace center bead to show simple harmonic motion
  int mid = NB / 2;
  pil_add_trace(w, pil_on(b[mid], 0, 0), 0.05, 200);

  pil_element_hide(
      pil_add_spring(w, pil_at(0, y0), pil_on(b[0], 0, 0), K, 0.0));
  for (int j = 0; j + 1 < NB; j++) {
    pil_element_hide(
        pil_add_spring(w, pil_on(b[j], 0, 0), pil_on(b[j + 1], 0, 0), K, 0.0));
  }

  pil_element_hide(
      pil_add_spring(w, pil_on(b[NB - 1], 0, 0), pil_at(LEN, y0), K, 0.0));
  pil_add_ground(w, (pil_vec2){0, y0 + 0.5}, (pil_vec2){0, y0 - 0.5});
  pil_add_ground(w, (pil_vec2){LEN, y0 + 0.5}, (pil_vec2){LEN, y0 - 0.5});
}

int main(int argc, char **argv) {
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);

  chain(w, 2.4, 1);
  chain(w, 0.0, 2);
  chain(w, -2.4, 3);
  pil_world_add_overlay(w, draw_chains, NULL);

  pil_add_label(w, pil_at(-0.4, 3.4), "y_j propto sin(n pi j / (N+1))", 0,
                0);
  pil_add_label(w, pil_at(-0.4, -3.4),
                "omega_n = 2 sqrt{k/m} sin(n pi / 2(N+1))", 0, 0);

  for (int n = 1; n <= 3; n++) {
    char buf[64];
    double om = 2.0 * sqrt(K / M) * sin(n * M_PI / (2.0 * (NB + 1)));
    snprintf(buf, sizeof buf, "n = %d,  \xcf\x89 = %.2f", n, om);
    pil_add_label(w, pil_at(LEN + 0.2, 2.4 - 2.4 * (n - 1)), buf, 2, 0);

    // Add node/antinode text annotations for clarity
    if (n == 1) {
      pil_add_label(w, pil_at(LEN / 2.0, 2.4 + 0.6), "antinode", 0.5, 0);
    } else if (n == 2) {
      pil_add_label(w, pil_at(LEN / 2.0, 0.0 + 0.6), "node", 0.5, 0);
    } else if (n == 3) {
      pil_add_label(w, pil_at(LEN / 6.0, -2.4 - 0.6), "node", 0.5, 0);
      pil_add_label(w, pil_at(LEN / 3.0, -2.4 - 0.6), "antinode", 0.5, 0);
    }
  }

  pil_view view;
  pil_view_fit_region(&view, -0.4, -3.6, 12.2, 3.6, 1500, 800, 0.0, 0.0, 1.0,
                      1.0, 12);
  int rc = example_main(w, &view, "Standing waves on a bead chain", argc, argv,
                        "waves.png");

  pil_world_free(w);

  return rc;
}
