/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Harshit-Dhanwalkar
 */

/* pil_optics.h - geometrical optics at a plane interface, prisms, thin lenses
 */
#ifndef pil_OPTICS_H
#define pil_OPTICS_H

#ifdef __cplusplus
extern "C" {
#endif

/* Snell's law n1 sin(ti) = n2 sin(tt), angles from normal in radians
 * Returns 1 and sets *tt, or 0 on total internal reflection */
int pil_snell(double n1, double n2, double ti, double *tt);
/* Critical angle asin(n2 / n1) for n1 > n2 (NaN otherwise) and Brewster's
 * angle atan(n2 / n1) */
double pil_critical_angle(double n1, double n2);
double pil_brewster_angle(double n1, double n2);

/* Refract a unit direction d travelling in medium n1 into medium n2 across an
 * interface whose unit normal nrm points back towards n1 side (so d . nrm < 0)
 * Writes refracted unit direction to out and returns 1, or reflected direction
 * and returns 0 on total internal reflection */
int pil_refract(const double d[2], const double nrm[2], double n1, double n2,
                double out[2]);
void pil_reflect(const double d[2], const double nrm[2], double out[2]);

/* Fresnel power reflectances for light polarised perpendicular (s) and
 * parallel (p) to plane of incidence; 1 beyond critical angle */
double pil_fresnel_rs(double n1, double n2, double ti);
double pil_fresnel_rp(double n1, double n2, double ti);

/* Prism of apex angle `apex` and index n in air: deviation of a ray with
 * incidence angle ti on first face (NaN if it is totally reflected) and minimum
 * deviation 2 asin(n sin(apex / 2)) - apex */
double pil_prism_deviation(double apex, double n, double ti);
double pil_prism_min_deviation(double apex, double n);

/* Thin lens, 1/f = 1/do + 1/di (real object distance do > 0): image distance
 * di (negative = virtual, on object side) and lateral magnification -di / do
 * do = f gives di = infinity */
double pil_lens_image(double f, double dobj);
double pil_lens_magnification(double f, double dobj);

/* Ray through a thin lens at height y above axis with slope dy/dx = s:
 * slope after lens is s - y / f (a ray through centre is undeflected, a
 * ray parallel to axis goes through focus) */
double pil_lens_ray_slope(double f, double y, double s);

/* 3 principal rays from an object tip at (-dobj, h): returns slope before lens
 * of (0) ray parallel to axis, (1) ray through lens centre, (2) ray aimed
 * through (or, for dobj < f, away from) front focus; all start at tip */
void pil_lens_principal_slopes(double f, double dobj, double h, double s[3]);

/* Cauchy dispersion n(lambda) = a + b / lambda^2 (lambda in micrometres) and an
 * approximate display colour (rgb in 0..1) of a wavelength in nanometres
 * (380..780) */
double pil_cauchy_index(double a, double b, double lambda_um);
void pil_wavelength_rgb(double lambda_nm, double rgb[3]);

/* Trace a ray through a convex glass polygon of index n in air
 * Polygon has nv vertices v[i] in counter-clockwise order
 * Ray starts at `start` (outside) moving along `dir`
 * At every boundary ray refracts out, or reflects totally if it cannot leave
 * Polyline (start, boundary hits, ...) goes to pts[0..max_pts-1] with last
 * point `tail` units beyond last exit along final direction; returns number of
 * points (0 if ray misses glass) and final direction in *out_dir (may be NULL)
 */
int pil_glass_trace(const double (*v)[2], int nv, double n,
                    const double start[2], const double dir[2], double tail,
                    double (*pts)[2], int max_pts, double out_dir[2]);

#ifdef __cplusplus
}
#endif

#endif
