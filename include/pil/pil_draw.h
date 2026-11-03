/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/*
 * pil_draw.h - Cairo renderer that draws a pil_world in style of a textbook
 * figure (grey bodies, coil springs, dashpots, hatched ground, hinges,
 * dimension lines, italic serif labels)
 */

#ifndef PIL_DRAW_H
#define PIL_DRAW_H

#include "pil_core.h"
#include <cairo.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Mapping world (metres, y up) <-> device pixels/points. */
struct pil_view {
  double cx, cy;        /* world point shown at centre of picture */
  double scale;         /* device units per metre                 */
  double width, height; /* device size                            */
};

typedef struct {
  double line_w;         /* body outline width                      */
  double thin_w;         /* springs, dimensions, thin lines         */
  double body_fill;      /* default grey (0..1)                     */
  double font_px;        /* label size                              */
  const char *font;      /* text family (italic)                    */
  const char *math_font; /* family for Greek letters and script l   */
  double ink[3];         /* line colour                             */
  double paper[3];       /* background                              */
  int draw_background;
  int text_halo; /* paper-coloured outline behind labels (readable over lines)
                  */
  int show_com;  /* draw centre-of-mass symbol                     */
  int spring_coils;   /* default number of coils                        */
  int upright;        /* set: Latin text in upright face (plot numbers) */
  double spring_w_px; /* default coil half-width in px                  */
  double hinge_r_px;
} pil_style;

pil_style pil_style_textbook(void);

void pil_view_fit(pil_view *v, double xmin, double ymin, double xmax,
                  double ymax, double width, double height, double margin_px);
/* Like pil_view_fit, but fits box into a sub-rectangle of canvas (fractions,
 * origin top-left) */
void pil_view_fit_region(pil_view *v, double xmin, double ymin, double xmax,
                         double ymax, double width, double height, double fx0,
                         double fy0, double fx1, double fy1, double margin_px);
pil_vec2 pil_view_to_screen(const pil_view *v, pil_vec2 world);
pil_vec2 pil_view_to_world(const pil_view *v, pil_vec2 screen);

/* Draw everything in world */
void pil_draw_world(cairo_t *cr, const pil_world *w, const pil_view *v,
                    const pil_style *st);

/* Export current state of world
 * Return 0 on success */
int pil_export_png(const pil_world *w, const pil_view *v, const pil_style *st,
                   const char *path);
int pil_export_pdf(const pil_world *w, const pil_view *v, const pil_style *st,
                   const char *path);
int pil_export_svg(const pil_world *w, const pil_view *v, const pil_style *st,
                   const char *path);

/* Simulate and record to a video/GIF through ffmpeg (must be installed)
 * Extension decides container (.mp4, .gif, .webm)
 * Returns 0 on success */
int pil_record(pil_world *w, const pil_view *v, const pil_style *st,
               const char *path, double seconds, double fps, double max_dt);

/* Same, with slow motion: one second of video shows `time_scale` seconds of
 * simulated time (time_scale = 0.001 animates a millisecond-scale circuit) */
int pil_record_ex(pil_world *w, const pil_view *v, const pil_style *st,
                  const char *path, double seconds, double fps, double max_dt,
                  double time_scale);

/* public primitives (screen coordinates, y down) */
void pil_draw_spring(cairo_t *cr, pil_vec2 a, pil_vec2 b, int coils,
                     double half_width, const pil_style *st);
void pil_draw_dashpot(cairo_t *cr, pil_vec2 a, pil_vec2 b,
                      double nominal_len_px, double half_width,
                      const pil_style *st);
void pil_draw_hinge(cairo_t *cr, pil_vec2 p, const pil_style *st);
void pil_draw_com(cairo_t *cr, pil_vec2 p, const pil_style *st);
void pil_draw_hatch(cairo_t *cr, pil_vec2 a, pil_vec2 b,
                    const pil_style *st); /* hatching on right of a->b */
void pil_draw_arrow(cairo_t *cr, pil_vec2 from, pil_vec2 to, double head_px,
                    const pil_style *st);
/* Text with x_1, x_{12}, x^2 sub/superscripts. align: 0 left, 0.5 centre, 1
 * right */
void pil_draw_text(cairo_t *cr, pil_vec2 p, const char *text, double align,
                   const pil_style *st);

#ifdef __cplusplus
}
#endif
#endif
