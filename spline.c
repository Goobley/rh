/* ------- file: -------------------------- spline.c ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Feb 16 14:57:56 1999   --

       --------------------------                      ----------RH-- */

/* --- Cubic spline interpolation routines. --         -------------- */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "rh.h"

void init_SplineState(SplineState* s)
{
  memset(s, 0, sizeof(SplineState));
  // *s = {NULL};
}

void free_SplineState(SplineState* s)
{
  if (s->M)
  {
    free(s->M);
    s->M = NULL;
  }
  if (s->u)
  {
    free(s->u);
    s->u = NULL;
  }
}

void spline_coef_stateless(SplineState* s, int N, double* x, double* y)
{
  s->ascend = (x[1] > x[0]) ? TRUE : FALSE;
  s->xmin = (s->ascend) ? x[0] : x[N - 1];
  s->xmax = (s->ascend) ? x[N - 1] : x[0];

  double* q = s->M = (double*)realloc(s->M, N * sizeof(double));
  s->u = (double*)realloc(s->u, N * sizeof(double));
  double hj = x[1] - x[0];
  double D = (y[1] - y[0]) / hj;

  q[0] = s->u[0] = 0.0;
  for (int j = 1; j < N - 1; j++) {
    double hj1 = x[j + 1] - x[j];
    double mu = hj / (hj + hj1);
    double D1 = (y[j + 1] - y[j]) / hj1;

    double p = mu * q[j - 1] + 2;
    q[j] = (mu - 1) / p;
    s->u[j] = ((D1 - D) * 6 / (hj + hj1) - mu * s->u[j - 1]) / p;

    hj = hj1;
    D = D1;
  }

  s->M[N - 1] = 0.0;
  for (int j = N - 2; j >= 0; j--) {
    s->M[j] = q[j] * s->M[j + 1] + s->u[j];
  }
  s->Ntable = N;
  s->xtable = x;
  s->ytable = y;
}

void spline_eval_stateless(SplineState* s, int N, double *x, double *y, bool_t hunt)
{

  for (int n = 0; n < N; n++) {
    if (x[n] <= s->xmin)
      y[n] = (s->ascend) ? s->ytable[0] : s->ytable[s->Ntable - 1];
    else if (x[n] >= s->xmax)
      y[n] = (s->ascend) ? s->ytable[s->Ntable - 1] : s->ytable[0];
    else {
      int j = 0;
      if (hunt)
        Hunt(s->Ntable, s->xtable, x[n], &j);
      else
        Locate(s->Ntable, s->xtable, x[n], &j);

      double hj = s->xtable[j + 1] - s->xtable[j];
      double fx = (x[n] - s->xtable[j]) / hj;
      double fx1 = 1 - fx;

      y[n] = fx1 * s->ytable[j] + fx * s->ytable[j + 1] +
             (fx1 * (SQ(fx1) - 1) * s->M[j] + fx * (SQ(fx) - 1) * s->M[j + 1]) *
                 SQ(hj) / 6.0;
    }
  }
}
/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

static bool_t ascend;
static int Ntable;
static double *xtable, xmin, xmax;

/* ------- begin -------------------------- splineCoef.c ------------ */

static double *M = NULL, *ytable;

void splineCoef(int N, double *x, double *y) {
  register int j;
  static double *u = NULL;

  double p, *q, hj, hj1, D, D1, mu;

  ascend = (x[1] > x[0]) ? TRUE : FALSE;
  xmin = (ascend) ? x[0] : x[N - 1];
  xmax = (ascend) ? x[N - 1] : x[0];

  q = M = (double *)realloc(M, N * sizeof(double));
  u = (double *)realloc(u, N * sizeof(double));
  hj = x[1] - x[0];
  D = (y[1] - y[0]) / hj;

  q[0] = u[0] = 0.0;
  for (j = 1; j < N - 1; j++) {
    hj1 = x[j + 1] - x[j];
    mu = hj / (hj + hj1);
    D1 = (y[j + 1] - y[j]) / hj1;

    p = mu * q[j - 1] + 2;
    q[j] = (mu - 1) / p;
    u[j] = ((D1 - D) * 6 / (hj + hj1) - mu * u[j - 1]) / p;

    hj = hj1;
    D = D1;
  }

  M[N - 1] = 0.0;
  for (j = N - 2; j >= 0; j--) {
    M[j] = q[j] * M[j + 1] + u[j];
  }
  Ntable = N;
  xtable = x;
  ytable = y;
}
/* ------- end ---------------------------- splineCoef.c ------------ */

/* ---------------------------------------- splineEval.c ------------ */

void splineEval(int N, double *x, double *y, bool_t hunt) {
  register int n;

  int j = 0;
  double hj, fx, fx1;

  for (n = 0; n < N; n++) {
    if (x[n] <= xmin)
      y[n] = (ascend) ? ytable[0] : ytable[Ntable - 1];
    else if (x[n] >= xmax)
      y[n] = (ascend) ? ytable[Ntable - 1] : ytable[0];
    else {
      if (hunt)
        Hunt(Ntable, xtable, x[n], &j);
      else
        Locate(Ntable, xtable, x[n], &j);

      hj = xtable[j + 1] - xtable[j];
      fx = (x[n] - xtable[j]) / hj;
      fx1 = 1 - fx;

      y[n] = fx1 * ytable[j] + fx * ytable[j + 1] +
             (fx1 * (SQ(fx1) - 1) * M[j] + fx * (SQ(fx) - 1) * M[j + 1]) *
                 SQ(hj) / 6.0;
    }
  }
}
/* ------- end ---------------------------- splineEval.c ------------ */
