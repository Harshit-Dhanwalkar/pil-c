/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* app_sdl.c - interactive window: SDL2 shows what Cairo draws */

#include "pil/pil_app.h"
#include "pil_internal.h"
#include <SDL.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

pil_app_config pil_app_defaults(const char *title, int width, int height) {
  pil_app_config c;
  memset(&c, 0, sizeof c);

  c.title = title;
  c.width = width;
  c.height = height;
  c.style = pil_style_textbook();
  c.view.cx = 0;
  c.view.cy = 0;
  c.view.scale = 100;
  c.view.width = width;
  c.view.height = height;
  c.time_scale = 1.0;
  c.max_dt = 5e-4;
  c.grab_k = 300.0;
  c.show_hud = 1;

  return c;
}

static void draw_hud(cairo_t *cr, const pil_world *w, const pil_app_config *c,
                     int paused, double speed) {
  char line[256];
  snprintf(line, sizeof line,
           "t = %.2f s   E = %.3f J   |err| = %.1e   x%.2f %s     [space] "
           "pause  [r] reset  [+/-] speed  [p/f/g] png/pdf/svg  [h] hide",
           pil_world_time(w), pil_world_energy(w),
           pil_world_constraint_error(w), speed, paused ? " PAUSED" : "");
  cairo_save(cr);
  cairo_select_font_face(cr, "monospace", CAIRO_FONT_SLANT_NORMAL,
                         CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, 12);
  cairo_set_source_rgba(cr, 0.35, 0.35, 0.35, 0.9);
  cairo_move_to(cr, 8, c->view.height - 8);
  cairo_show_text(cr, line);
  cairo_restore(cr);
}

int pil_app_run(pil_world *w, pil_app_config *cfg) {
  if (SDL_Init(SDL_INIT_VIDEO) != 0) {
    fprintf(stderr, "SDL_Init: %s\n", SDL_GetError());
    return -1;
  }
  SDL_Window *win = SDL_CreateWindow(
      cfg->title ? cfg->title : "pil", SDL_WINDOWPOS_CENTERED,
      SDL_WINDOWPOS_CENTERED, cfg->width, cfg->height, SDL_WINDOW_RESIZABLE);
  SDL_Renderer *ren =
      win ? SDL_CreateRenderer(win, -1, SDL_RENDERER_PRESENTVSYNC) : NULL;
  if (!ren && win)
    ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_SOFTWARE);
  if (!win || !ren) {
    fprintf(stderr, "SDL window: %s\n", SDL_GetError());
    SDL_Quit();
    return -1;
  }

  int W = cfg->width, H = cfg->height;
  SDL_Texture *tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_ARGB8888,
                                       SDL_TEXTUREACCESS_STREAMING, W, H);
  cairo_surface_t *surf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, W, H);
  cairo_t *cr = cairo_create(surf);

  pil_world_snapshot(w);
  int paused = cfg->start_paused, running = 1, frames = 0, panning = 0,
      step_once = 0;
  double speed = cfg->time_scale > 0 ? cfg->time_scale : 1.0;
  Uint64 last = SDL_GetPerformanceCounter();
  const double freq = (double)SDL_GetPerformanceFrequency();
  int px = 0, py = 0;

  while (running) {
    SDL_Event e;
    while (SDL_PollEvent(&e)) {
      switch (e.type) {
      case SDL_QUIT:
        running = 0;
        break;
      case SDL_WINDOWEVENT:
        if (e.window.event == SDL_WINDOWEVENT_SIZE_CHANGED) {
          W = e.window.data1;
          H = e.window.data2;
          if (W < 16 || H < 16) {
            break;
          }
          cairo_destroy(cr);
          cairo_surface_destroy(surf);
          SDL_DestroyTexture(tex);
          tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_ARGB8888,
                                  SDL_TEXTUREACCESS_STREAMING, W, H);
          surf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, W, H);
          cr = cairo_create(surf);
          cfg->view.width = W;
          cfg->view.height = H;
        }
        break;
      case SDL_KEYDOWN:
        switch (e.key.keysym.sym) {
        case SDLK_ESCAPE:
        case SDLK_q:
          running = 0;
          break;
        case SDLK_SPACE:
          paused = !paused;
          break;
        case SDLK_r:
          pil_world_reset(w);
          break;
        case SDLK_RIGHT:
          step_once = 1;
          break;
        case SDLK_PLUS:
        case SDLK_EQUALS:
        case SDLK_KP_PLUS:
          speed *= 1.25;
          break;
        case SDLK_MINUS:
        case SDLK_KP_MINUS:
          speed /= 1.25;
          break;
        case SDLK_h:
          cfg->show_hud = !cfg->show_hud;
          break;
        case SDLK_p:
          pil_export_png(w, &cfg->view, &cfg->style, "pil_frame.png");
          fprintf(stderr, "saved pil_frame.png\n");
          break;
        case SDLK_f:
          pil_export_pdf(w, &cfg->view, &cfg->style, "pil_frame.pdf");
          fprintf(stderr, "saved pil_frame.pdf\n");
          break;
        case SDLK_g:
          pil_export_svg(w, &cfg->view, &cfg->style, "pil_frame.svg");
          fprintf(stderr, "saved pil_frame.svg\n");
          break;
        default:
          break;
        }
        break;
      case SDL_MOUSEBUTTONDOWN:
        if (e.button.button == SDL_BUTTON_LEFT) {
          pil_vec2 wp = pil_view_to_world(&cfg->view,
                                          pil_v2(e.button.x, e.button.y)),
                   loc;
          pil_body *b = pil_world_pick(w, wp, &loc);
          if (b)
            pil_world_set_grab(w, b, loc, wp, cfg->grab_k * b->mass,
                               2.0 * sqrt(cfg->grab_k) * b->mass);
        } else {
          panning = 1;
          px = e.button.x;
          py = e.button.y;
        }
        break;
      case SDL_MOUSEBUTTONUP:
        if (e.button.button == SDL_BUTTON_LEFT)
          pil_world_clear_grab(w);
        else
          panning = 0;
        break;
      case SDL_MOUSEMOTION:
        if (w->grab.b)
          w->grab.target =
              pil_view_to_world(&cfg->view, pil_v2(e.motion.x, e.motion.y));
        if (panning) {
          cfg->view.cx -= (e.motion.x - px) / cfg->view.scale;
          cfg->view.cy += (e.motion.y - py) / cfg->view.scale;
          px = e.motion.x;
          py = e.motion.y;
        }
        break;
      case SDL_MOUSEWHEEL: {
        int mx;
        int my;
        SDL_GetMouseState(&mx, &my);
        pil_vec2 before = pil_view_to_world(&cfg->view, pil_v2(mx, my));
        cfg->view.scale *= pow(1.1, e.wheel.y);
        pil_vec2 after = pil_view_to_world(&cfg->view, pil_v2(mx, my));
        cfg->view.cx += before.x - after.x;
        cfg->view.cy += before.y - after.y;
        break;
      }
      default:
        break;
      }
    }

    Uint64 now = SDL_GetPerformanceCounter();
    double dt = (double)(now - last) / freq;
    last = now;
    if (dt > 0.05) {
      dt = 0.05; // never "catch up" after a stall
    }
    if (!paused) {
      pil_world_advance(w, dt * speed, cfg->max_dt);
    } else if (step_once) {
      pil_world_advance(w, 1.0 / 60.0, cfg->max_dt);
    }
    step_once = 0;

    if (cfg->on_frame) {
      cfg->on_frame(w, pil_world_time(w), cfg->user);
    }
    pil_draw_world(cr, w, &cfg->view, &cfg->style);
    if (cfg->show_hud) {
      draw_hud(cr, w, cfg, paused, speed);
    }
    cairo_surface_flush(surf);
    SDL_UpdateTexture(tex, NULL, cairo_image_surface_get_data(surf),
                      cairo_image_surface_get_stride(surf));
    SDL_RenderClear(ren);
    SDL_RenderCopy(ren, tex, NULL, NULL);
    SDL_RenderPresent(ren);

    if (cfg->max_frames && ++frames >= cfg->max_frames) {
      running = 0;
    }
  }

  if (cfg->exit_screenshot) {
    pil_export_png(w, &cfg->view, &cfg->style, cfg->exit_screenshot);
  }

  cairo_destroy(cr);
  cairo_surface_destroy(surf);
  SDL_DestroyTexture(tex);
  SDL_DestroyRenderer(ren);
  SDL_DestroyWindow(win);
  SDL_Quit();

  return 0;
}
