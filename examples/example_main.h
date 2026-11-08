/* Shared entry point for examples:
 *   prog [out.png]              render initial frame
 *   prog --record out.gif [s]   simulate s seconds (default 4) into a gif/mp4
 *   prog --live                 interactive window (needs SDL2)
 */

#ifndef EXAMPLE_MAIN_H
#define EXAMPLE_MAIN_H
#include <pil/pil_core.h>
#include <pil/pil_draw.h>
#ifdef PIL_HAVE_SDL
#include <pil/pil_app.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* An example can slow itself down (circuits live in milliseconds): define
 * PIL_EXAMPLE_TIME_SCALE before including this header */
// Sim time to run before the still image is drawn, for figures that only show
// their currents and voltages once running (circuits)
#ifndef PIL_EXAMPLE_IMAGE_TIME
#define PIL_EXAMPLE_IMAGE_TIME 0.0
#endif

#ifndef PIL_EXAMPLE_TIME_SCALE
#define PIL_EXAMPLE_TIME_SCALE 1.0
#endif

static int example_main(pil_world *w, pil_view *view, const char *title,
                        int argc, char **argv, const char *default_png) {
  pil_style st = pil_style_textbook();
  int rc = 0;
  (void)title; // only used by SDL window

  // dots are specified in px per second of wall time; sim runs slower or faster
  // by time scale
  st.current_dot_speed /= PIL_EXAMPLE_TIME_SCALE;

  if (argc > 1 && !strcmp(argv[1], "--live")) {
#ifdef PIL_HAVE_SDL
    pil_app_config cfg =
        pil_app_defaults(title, (int)view->width, (int)view->height);
    cfg.view = *view;
    cfg.style = st;
    cfg.time_scale = PIL_EXAMPLE_TIME_SCALE;

    for (int i = 2; i + 1 < argc; i++) {
      if (!strcmp(argv[i], "--frames")) {
        cfg.max_frames = atoi(argv[i + 1]);
      }
      if (!strcmp(argv[i], "--shot")) {
        cfg.exit_screenshot = argv[i + 1];
      }
    }

    rc = pil_app_run(w, &cfg);
#else
    fprintf(stderr, "built without SDL2: --live unavailable\n");
    rc = 1;
#endif
  } else if (argc > 2 && !strcmp(argv[1], "--record")) {
    double secs = argc > 3 ? atof(argv[3]) : 4.0;
    rc = pil_record_ex(w, view, &st, argv[2], secs, 30.0, 5e-4,
                       PIL_EXAMPLE_TIME_SCALE);
    printf(rc ? "record failed (is ffmpeg installed?)\n"
              : "recorded %s (%.1f s)\n",
           argv[2], secs);
  } else if (argc > 1 && !strcmp(argv[1], "--check")) {
    double secs = argc > 2 ? atof(argv[2]) : 10.0;
    double E0 = pil_world_energy(w);
    double maxerr = 0;
    for (double t = 0; t < secs; t += 0.01) {
      pil_world_advance(w, 0.01, 5e-4);
      double e = pil_world_constraint_error(w);
      if (e > maxerr) {
        maxerr = e;
      }
    }

    double E1 = pil_world_energy(w);
    printf("t = %.1f s   E: %.6g -> %.6g   max |constraint error| = %.2e\n",
           secs, E0, E1, maxerr);
    rc = (E1 == E1 && maxerr < 1e-3) ? 0 : 1; // NaN or a large error fails
    if (rc) {
      printf("CHECK FAILED\n");
    }
  } else {
    const char *png = argc > 1 ? argv[1] : default_png;
    if (PIL_EXAMPLE_IMAGE_TIME > 0) {
      pil_world_advance(w, PIL_EXAMPLE_IMAGE_TIME, 5e-4);
    }

    rc = pil_export_png(w, view, &st, png);
    printf(rc ? "export failed\n" : "wrote %s\n", png);
  }

  return rc;
}
#endif
