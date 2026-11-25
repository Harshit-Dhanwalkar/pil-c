#ifndef PIL_TIMER555_H
#define PIL_TIMER555_H

#ifdef __cplusplus
extern "C" {
#endif

/* 555 timer in astable mode: C charges through R1 + R2 towards Vcc until
 * v_C = 2/3 Vcc, then discharges through R2 until v_C = 1/3 Vcc; output is
 * high while charging */
typedef struct {
  double r1;
  double r2;
  double c;
  double vcc;
  double vc; /* Capacitor voltage */
  int out;   /* 1 = output high (charging) */
} pil_555;

void pil_555_init(pil_555 *t, double r1, double r2, double c, double vcc);

/* advance by dt using exact exponential solution, with threshold crossing
 * resolved inside step */
void pil_555_step(pil_555 *t, double dt);

double pil_555_high_time(double r1, double r2, double c); /* ln2 (R1+R2) C */
double pil_555_low_time(double r1, double r2, double c);  /* ln2 R2 C */
double pil_555_frequency(double r1, double r2, double c);
double pil_555_duty(double r1, double r2); /* (R1+R2)/(R1+2R2) */

#ifdef __cplusplus
}
#endif

#endif
