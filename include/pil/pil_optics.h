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

#ifdef __cplusplus
}
#endif

#endif
