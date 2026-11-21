/*
 * doppler.c - Doppler effect of a moving source
 * A source moving at v = c / 2 emits a pulse every T_0; each pulse spreads as
 * a circle about place where it was emitted. Fronts bunch up ahead (wave
 * length (1 - v/c) c T_0) and spread out behind ((1 + v/c) c T_0)
 * Plot counts pulses received by an observer ahead and one behind: slopes are
 * f / (1 - v/c) and f / (1 + v/c), a ratio of 3 here
 */

#define PIL_EXAMPLE_TIME_SCALE 1.0
#define PIL_EXAMPLE_IMAGE_TIME 12.0
#include "example_main.h"
#include <math.h>

#define C_ 1.0
#define V_ 0.5
#define T0_ 0.6
#define X0_ -3.0
#define XA_ 5.5  // Observer ahead
#define XB_ -5.5 // Observer behind

// Pulses received by an observer at xo by time t
static int received(double xo, double t) {
  int n = 0;
  for (int i = 0; i * T0_ <= t; i++) {
    double te = i * T0_;
    double xs = X0_ + V_ * te;
    if (te + fabs(xo - xs) / C_ <= t) {
      n++;
    }
  }

  return n;
}

static double sig_ahead(const pil_world *w, void *u) {
  (void)u;

  return received(XA_, pil_world_time(w));
}

static double sig_behind(const pil_world *w, void *u) {
  (void)u;

  return received(XB_, pil_world_time(w));
}

static void overlay(cairo_t *cr, const pil_view *v, void *user) {
  const pil_world *w = user;
  pil_style st = pil_style_textbook();
  double t = pil_world_time(w);

  double sc = pil_view_to_screen(v, pil_v2(1, 0)).x -
              pil_view_to_screen(v, pil_v2(0, 0)).x;
  double cx1;
  double cy1;
  double cx2;
  double cy2;
  cairo_clip_extents(cr, &cx1, &cy1, &cx2, &cy2);
  cairo_save(cr);
  cairo_rectangle(cr, cx1, cy1, cx2 - cx1,
                  0.66 * (cy2 - cy1)); // above the plot
  cairo_clip(cr);
  cairo_set_line_width(cr, st.thin_w * 1.2);
  for (int i = 0; i * T0_ <= t; i++) {
    double te = i * T0_; double r = C_ * (t - te);
    if (r > 12.0) {
      continue;
    }

    pil_vec2 o = pil_view_to_screen(v, pil_v2(X0_ + V_ * te, 0));
    double g = 0.1 + 0.6 * r / 12.0;
    cairo_set_source_rgb(cr, g, g, g);
    cairo_arc(cr, o.x, o.y, r * sc, 0, 2 * M_PI);
    cairo_stroke(cr);
  }

  cairo_restore(cr);

  // Source and its velocity arrow
  pil_vec2 s = pil_view_to_screen(v, pil_v2(X0_ + V_ * t, 0));
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_arc(cr, s.x, s.y, 7, 0, 2 * M_PI);
  cairo_fill(cr);
  cairo_set_line_width(cr, st.line_w);
  cairo_move_to(cr, s.x + 10, s.y);
  cairo_line_to(cr, s.x + 10 + 0.8 * sc, s.y);
  cairo_stroke(cr);
  cairo_move_to(cr, s.x + 10 + 0.8 * sc, s.y);
  cairo_line_to(cr, s.x + 10 + 0.8 * sc - 8, s.y - 4);
  cairo_line_to(cr, s.x + 10 + 0.8 * sc - 8, s.y + 4);
  cairo_close_path(cr);
  cairo_fill(cr);
  pil_draw_text(cr, pil_v2(s.x + 14 + 0.4 * sc, s.y - 14), "v", 0.5, &st);

  // Observers
  double xs[2] = {XA_, XB_};
  const char *nm[2] = {"A", "B"};
  for (int k = 0; k < 2; k++) {
    pil_vec2 q = pil_view_to_screen(v, pil_v2(xs[k], 0));
    cairo_rectangle(cr, q.x - 6, q.y - 6, 12, 12);
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_fill_preserve(cr);
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_stroke(cr);
    pil_draw_text(cr, pil_v2(q.x, q.y + 24), nm[k], 0.5, &st);
  }
}

int main(int argc, char **argv) {
  pil_world *w = pil_world_new();
  pil_world_set_gravity(w, 0, 0);
  pil_world_add_overlay(w, overlay, w);

  pil_plot *p =
      pil_add_plot(w, 0.0, 0.68, 0.62, 0.99, "t [s]", "pulses received", 14.0);
  pil_plot_add_signal(p, "A (ahead)", sig_ahead, NULL);
  pil_plot_add_signal(p, "B (behind)", sig_behind, NULL);
  pil_plots_record(w, 0.02, 1e-2);

  pil_view view;
  pil_view_fit_region(&view, -7.0, -4.0, 7.0, 4.0, 1300, 760, 0.0, 0.0, 1.0,
                      0.66, 14);
  int rc = example_main(w, &view, "Doppler effect", argc, argv, "doppler.png");

  pil_world_free(w);

  return rc;
}
