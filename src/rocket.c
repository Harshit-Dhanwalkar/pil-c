#include <math.h>
#include <pil/pil_ode.h>
#include <pil/pil_rocket.h>

void pil_rocket_init(pil_rocket_state *s) {
  s->t = 0.0;
  s->h = 0.0;
  s->v = 0.0;
  s->stage = 0;
  s->mf = -1.0;
}

static void ensure_fuel(const pil_rocket *r, pil_rocket_state *s) {
  if (s->mf < 0.0) {
    s->mf = s->stage < r->n ? r->st[s->stage].m_fuel : 0.0;
  }
}

double pil_rocket_mass(const pil_rocket *r, const pil_rocket_state *s) {
  double m = r->payload;
  for (int i = s->stage; i < r->n; i++) {
    m += r->st[i].m_dry +
         (i == s->stage && s->mf >= 0.0 ? s->mf : r->st[i].m_fuel);
  }

  return m;
}

int pil_rocket_burning(const pil_rocket *r, const pil_rocket_state *s) {
  return s->stage < r->n;
}

double pil_rocket_dv(double ve, double m0, double m1) {
  return ve * log(m0 / m1);
}

typedef struct {
  const pil_rocket *r;
  double m_fixed; /* Everything except propellant of firing stage */
  const pil_stage *st;
} burn_ctx;

static void burn_deriv(double t, const double *y, double *dydt, void *user) {
  (void)t;
  const burn_ctx *c = user;
  double m = c->m_fixed + y[2];
  dydt[0] = y[1];
  dydt[1] = c->st->mdot * c->st->ve / m - c->r->g;
  dydt[2] = -c->st->mdot;
}

static void coast_deriv(double t, const double *y, double *dydt, void *user) {
  (void)t;
  const pil_rocket *r = user;
  dydt[0] = y[1];
  dydt[1] = -r->g;
}

void pil_rocket_step(const pil_rocket *r, pil_rocket_state *s, double dt) {
  int guard = 0;
  while (dt > 0.0 && guard++ < 16) {
    if (s->stage >= r->n) {
      double y[2] = {s->h, s->v};
      pil_rk4_step(coast_deriv, 2, s->t, y, dt, (void *)r);
      s->h = y[0];
      s->v = y[1];
      s->t += dt;
      return;
    }

    ensure_fuel(r, s);
    const pil_stage *st = &r->st[s->stage];
    double h = dt, burn_time = s->mf / st->mdot;
    int burnout = 0;
    if (burn_time <= dt) {
      h = burn_time;
      burnout = 1;
    }

    burn_ctx c = {r, pil_rocket_mass(r, s) - s->mf, st};
    double y[3] = {s->h, s->v, s->mf};
    if (h > 0.0) {
      pil_rk4_step(burn_deriv, 3, s->t, y, h, &c);
    }
    s->h = y[0];
    s->v = y[1];
    s->mf = y[2];
    s->t += h;
    dt -= h;
    if (burnout) {
      s->stage++;
      s->mf = -1.0;
    }
  }
}
