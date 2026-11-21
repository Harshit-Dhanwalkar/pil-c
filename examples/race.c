/*
 * race.c - rolling down an incline: hoop, disc, sphere and a frictionless
 * slider Same mass and radius, same slope, same start. For rolling without
 * slipping acceleration along slope is
 *     a = g sin(alpha) / (1 + I / (m R^2))
 * so more of mass sits near rim (large I), slower body: hoop (I = m R^2), disc
 * (1/2), sphere (2/5); the frictionless slider (no rotation) has a = g
 * sin(alpha)
 * Plot shows the distance travelled along the slope
 * Everything is public-API
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 1.0
#include "example_main.h"
#include <math.h>

#define NLANE 4

typedef struct {
  pil_body *b;
  pil_vec2 start;
} racer;

static racer racers[NLANE];
static pil_vec2 slope_dir;

static double sig_s(const pil_world *w, void *u) {
  (void)w;
  const racer *r = u;
  pil_vec2 d = {r->b->pos.x - r->start.x, r->b->pos.y - r->start.y};

  return d.x * slope_dir.x + d.y * slope_dir.y;
}

int main(int argc, char **argv) {
  const double al = 25.0 * M_PI / 180.0, R = 0.2, len = 5.6, pitch = 0.75;
  const double kI[NLANE] = {1.0, 0.5, 0.4, -1.0}; // I / (m R^2); -1: slider
  const char *names[NLANE] = {"hoop  I = m R^2", "disc  I = 1/2 m R^2",
                              "sphere  I = 2/5 m R^2", "slider (frictionless)"};
  slope_dir = (pil_vec2){cos(al), -sin(al)};
  pil_vec2 nrm = {sin(al), cos(al)};
  pil_world *w = pil_world_new();

  for (int i = 0; i < NLANE; i++) {
    pil_vec2 p0 = {-1.6 + nrm.x * pitch * i, 3.3 + nrm.y * pitch * i - 0.6 * 0};
    p0.x -= 0.0;
    pil_vec2 p1 = {p0.x + slope_dir.x * len, p0.y + slope_dir.y * len};
    pil_add_ground(w, p0, p1);
    pil_vec2 c = {p0.x + slope_dir.x * 0.4 + nrm.x * R,
                  p0.y + slope_dir.y * 0.4 + nrm.y * R};
    pil_body_desc d = {.shape = pil_disc(R),
                       .mass = 1.0,
                       .pos = c,
                       .flags = PIL_SPOKE,
                       .fill = 0.85 - 0.18 * i};
    if (kI[i] > 0) {
      d.inertia = kI[i] * 1.0 * R * R;
    }

    racers[i].b = pil_add_body(w, &d);
    racers[i].start = c;
    if (kI[i] > 0) {
      pil_add_rolling(w, racers[i].b, (pil_vec2){0, 0}, R, p0, slope_dir);
    } else {
      pil_add_floor_contact(w, racers[i].b, pil_at(p0.x, p0.y), slope_dir, 3e5,
                            2e3, 0.0);
    }

    pil_add_label(w, pil_at(p1.x + 0.15, p1.y + 0.05), names[i], 2, 0);
  }

  pil_plot *pl = pil_add_plot(w, 0.08, 0.70, 0.98, 0.98, "t [s]", "s [m]", 1.6);
  for (int i = 0; i < NLANE; i++) {
    pil_plot_add_signal(
        pl,
        i == 3 ? "slider" : (i == 0 ? "hoop" : (i == 1 ? "disc" : "sphere")),
        sig_s, &racers[i]);
  }

  pil_plots_record(w, 0.02, 5e-4);
  pil_add_gravity_indicator(w, (pil_vec2){-1.9, 4.9});

  pil_view view;
  pil_view_fit_region(&view, -2.2, -0.9, 7.4, 5.4, 1500, 950, 0.0, 0.0, 1.0,
                      0.68, 12);
  int rc =
      example_main(w, &view, "Rolling down an incline", argc, argv, "race.png");

  pil_world_free(w);

  return rc;
}
