/* Shared entry point for the examples:
 *   prog [out.png]              render the initial frame
 *   prog --record out.gif [s]   simulate s seconds (default 4) into a gif/mp4
 *   prog --live                 interactive window (needs SDL2)
 */
#ifndef EXAMPLE_MAIN_H
#define EXAMPLE_MAIN_H
#include <pil/pil_core.h>
#include <pil/pil_draw.h>
#ifdef pil_HAVE_SDL
#include <pil/pil_app.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int example_main(pil_world *w, pil_view *view, const char *title, int argc, char **argv,
                        const char *default_png)
{
    pil_style st = pil_style_textbook();
    int rc = 0;
    if (argc > 1 && !strcmp(argv[1], "--live")) {
#ifdef pil_HAVE_SDL
        pil_app_config cfg = pil_app_defaults(title, (int)view->width, (int)view->height);
        cfg.view = *view; cfg.style = st;
        for (int i = 2; i + 1 < argc; i++) {
            if (!strcmp(argv[i], "--frames")) cfg.max_frames = atoi(argv[i + 1]);
            if (!strcmp(argv[i], "--shot")) cfg.exit_screenshot = argv[i + 1];
        }
        rc = pil_app_run(w, &cfg);
#else
        fprintf(stderr, "built without SDL2: --live unavailable\n"); rc = 1;
#endif
    } else if (argc > 2 && !strcmp(argv[1], "--record")) {
        double secs = argc > 3 ? atof(argv[3]) : 4.0;
        rc = pil_record(w, view, &st, argv[2], secs, 30.0, 5e-4);
        printf(rc ? "record failed (is ffmpeg installed?)\n" : "recorded %s (%.1f s)\n", argv[2], secs);
    } else {
        const char *png = argc > 1 ? argv[1] : default_png;
        rc = pil_export_png(w, view, &st, png);
        printf(rc ? "export failed\n" : "wrote %s\n", png);
    }
    return rc;
}
#endif
