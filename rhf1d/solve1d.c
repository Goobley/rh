/* ---------------------------------------- solve1d.c --------------- */
/* ------- One-D, plane-parallel version ---------------------------- */

/* Main routine of my_1D_radiative_transfer_program.
 *
 * Han Uitenbroek
 * Last modified: Wed Apr 22 09:44:52 2009 --
 */

#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "accelerate.h"
#include "constant.h"
#include "statistics.h"
#include "inputs.h"

/* --- Function prototypes --                          -------------- */

void GaussLeg(double x1, double x2, double *x, double *w, int n);

/* --- Global variables --                             -------------- */

enum Topology topology = ONE_D_PLANE;

int Nlambda;
double *Bp, *epsilon, *phi, *wlamb, wphi, BoundCond[4] = {0.0, 0.0, 1.0, 0.0},
                                          *chi, *Sny, **Iemerge;
char messageStr[MAX_LINE_SIZE];
struct Ng *NgS;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
CommandLine commandline;

/* ------- begin -------------------------- solve1d.c --------------- */

int main(int argc, char *argv[]) {
  register int l, la;

  int Ndep, Nrays, NmaxIter, Ngdelay, Ngorder, Ngperiod, Nread, Nwrite;
  double *lambda, Adamp, iterLimit;

  commandline.quiet = FALSE;
  commandline.logfile = stderr;

  stats.printCPU = TRUE;
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  /* --- Read input data --                        ------------------ */

  setbuf(stdout, NULL);

  geometry.vboundary[TOP] = ZERO;
  geometry.vboundary[BOTTOM] = REFLECTIVE;

  Nread = fread(&NmaxIter, sizeof(int), 1, stdin);
  Nread = fread(&iterLimit, sizeof(double), 1, stdin);
  Nread = fread(&Ngdelay, sizeof(int), 1, stdin);
  Nread = fread(&Ngorder, sizeof(int), 1, stdin);
  Nread = fread(&Ngperiod, sizeof(int), 1, stdin);

  Nread = fread(&geometry.Nrays, sizeof(int), 1, stdin);
  Nrays = geometry.Nrays;
  Nread = fread(&Nlambda, sizeof(int), 1, stdin);
  Nread = fread(&geometry.Ndep, sizeof(int), 1, stdin);
  Ndep = geometry.Ndep;

  geometry.height = (double *)malloc(Ndep * sizeof(double));
  Nread = fread(geometry.height, sizeof(double), Ndep, stdin);
  for (l = 0; l < Ndep; l++)
    geometry.height[l] *= KM_TO_M;

  lambda = (double *)malloc(Nlambda * sizeof(double));
  phi = (double *)malloc(Nlambda * sizeof(double));
  Nread = fread(&Adamp, sizeof(double), 1, stdin);
  Nread = fread(lambda, sizeof(double), Nlambda, stdin);
  wlamb = (double *)malloc(Nlambda * sizeof(double));
  wlamb[0] = (lambda[1] - lambda[0]);
  for (la = 1; la < Nlambda - 1; la++)
    wlamb[la] = (lambda[la + 1] - lambda[la - 1]);
  wlamb[Nlambda - 1] = (lambda[Nlambda - 1] - lambda[Nlambda - 2]);

  wphi = 0.0;
  for (la = 0; la < Nlambda; la++) {
    phi[la] = Voigt(Adamp, lambda[la], NULL, RYBICKI) / SQRTPI;
    wphi += wlamb[la] * phi[la];
  }
  wphi = 1.0 / wphi;

  chi = (double *)malloc(Ndep * sizeof(double));
  Bp = (double *)malloc(Ndep * sizeof(double));
  epsilon = (double *)malloc(Ndep * sizeof(double));
  Nread = fread(chi, sizeof(double), Ndep, stdin);
  Nread = fread(Bp, sizeof(double), Ndep, stdin);
  Nread = fread(epsilon, sizeof(double), Ndep, stdin);

  Sny = (double *)malloc(Ndep * sizeof(double));
  for (l = 0; l < Ndep; l++) {
    Sny[l] = Bp[l];
  }
  NgS = NgInit(Ndep, Ngdelay, Ngorder, Ngperiod, Sny);

  /* --- Determine angle values and weights --     ------------------ */

  geometry.wmu = (double *)malloc(Nrays * sizeof(double));
  geometry.muz = (double *)malloc(Nrays * sizeof(double));
  GaussLeg(0.0, 1.0, geometry.muz, geometry.wmu, Nrays);

  /* --- Iterate --                                ------------------ */

  Iterate(NmaxIter, iterLimit);

  /* --- Write output --                           ------------------ */

  Nwrite = fwrite(phi, sizeof(double), Nlambda, stdout);
  Nwrite = fwrite(&wphi, sizeof(double), 1, stdout);
  Nwrite = fwrite(Sny, sizeof(double), Ndep, stdout);
  Nwrite = fwrite(Iemerge[0], sizeof(double), Nrays * Nlambda, stdout);

  printTotalCPU();
}
/* ------- end ---------------------------- solve1d.c --------------- */
