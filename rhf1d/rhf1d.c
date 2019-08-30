/* ------- file: -------------------------- rhf1d.c -----------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Feb 24 16:40:14 2011 --

       --------------------------                      ----------RH-- */

/* --- Main routine of 1D plane-parallel radiative transfer program.
       MALI scheme formulated according to Rybicki & Hummer

  See: G. B. Rybicki and D. G. Hummer 1991, A&A 245, p. 171-181
       G. B. Rybicki and D. G. Hummer 1992, A&A 263, p. 209-215

       Formal solution is performed with Feautrier difference scheme
       in static atmospheres, and with piecewise quadratic integration
       in moving atmospheres.

       --                                              -------------- */

#include <string.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "statistics.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"

#ifdef NO_MAIN
  #define CMO_NO_PROF
#endif
#ifdef CMO_NO_PROF
  #undef CMO_NO_PROF
  #define SMALL_PROF
#endif
#define CMO_PROFILE_IMPL
#include "CmoProfile.h"


/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

#ifndef NO_MAIN
enum Topology topology = ONE_D_PLANE;

Atmosphere atmos = {0};
Geometry geometry = {0};
Spectrum spectrum = {0};
ProgramStats stats = {0};
InputData input = {0};
CommandLine commandline = {0};
char messageStr[MAX_MESSAGE_LENGTH];
#else
extern enum Topology topology;

extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern ProgramStats stats;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[];

#endif

static void init_profiler_id_sched(void* userdata, sched_uint threadId)
{
  CmoProfileId = threadId;
  // tss_create(&CmoProfileId, free);
  // tss_set(CmoProfileId, threadId)
}
/* ------- begin -------------------------- rhf1d.c ----------------- */
void cmo_init_misc()
{
  // return;
  printf("C Misc init\n");
  printf("C Misc init\n");
  printf("C Misc init\n");
  printf("C Misc init\n");
  printf("C Misc init\n");
  printf("C Misc init\n");
  printf("C Misc init\n");
  printf("C Misc init\n");
  fflush(stdout);
  // assert(FALSE);

  // FILE* f = fopen("blah.txt", "w");
  // fprintf(f, "test");
  // fflush(f);
  // fclose(f);

#ifndef SMALL_PROF
  cmo_allocate_prof_array(input.Nthreads, 128*1024*1024);
#else
  cmo_allocate_prof_array(input.Nthreads, 128);
#endif
  CMO_PROF_REGION_START("ProfileStart");
  CMO_PROF_REGION_END("ProfileStart");
  ProfileEntry* entry1 = &profiler[0].entries[profiler[0].end-2];
  ProfileEntry* entry2 = &profiler[0].entries[profiler[0].end-1];
  for (int i = 1; i < input.Nthreads; ++i)
  {
    cmo_prof_push_entry(&profiler[i], entry1->func, entry1->count, entry1->exit);
    cmo_prof_push_entry(&profiler[i], entry2->func, entry2->count, entry2->exit);
  }

  sched_size neededMem;
  scheduler_init(&input.sched, &neededMem, input.Nthreads, NULL);
  input.schedMem = calloc(neededMem, 1);
  scheduler_start(&input.sched, input.schedMem, init_profiler_id_sched);
  // cmo_prof_start_timed_region(&profiler[CmoProfileId], __func__);

  input.splineState = (SplineState*)calloc(input.Nthreads, sizeof(SplineState));
  for (int i = 0; i < input.Nthreads; ++i)
  {
    init_SplineState(&input.splineState[i]);
  }
}

void cmo_free_misc()
{
  CMO_PROF_REGION_START("ProfileEnd");
  CMO_PROF_REGION_END("ProfileEnd");
  ProfileEntry* entry1 = &profiler[0].entries[profiler[0].end-2];
  ProfileEntry* entry2 = &profiler[0].entries[profiler[0].end-1];
  for (int i = 1; i < input.Nthreads; ++i)
  {
    cmo_prof_push_entry(&profiler[i], entry1->func, entry1->count, entry1->exit);
    cmo_prof_push_entry(&profiler[i], entry2->func, entry2->count, entry2->exit);
  }
  for (int i = 0; i < input.Nthreads; ++i)
  {
    char buf[512];
    snprintf(buf, 512, "Prof_%02d.txt", i);
    cmo_prof_print_file(&profiler[i], buf);
  }
  // cmo_free_prof(&profiler[CmoProfileId]);
  cmo_free_prof_array(input.Nthreads);
  for (int i = 0; i < input.Nthreads; ++i)
  {
    free_SplineState(&input.splineState[i]);
  }
  free(input.splineState);
  input.splineState = NULL;
  if (input.Nthreads > 1)
  {
    scheduler_stop(&input.sched, 1);
    free(input.schedMem);
  }
  printTotalCPU();
}

#ifndef NO_MAIN
int main(int argc, char *argv[]) {
  bool_t write_analyze_output, equilibria_only;
  int niter, nact;

  Atom *atom;
  Molecule *molecule;

  /* --- Read input data and initialize --             -------------- */
  // cmo_init_prof(&threadProf, 0, 1024*32);
  // cmo_prof_start_timed_region(&threadProf, __func__);

  // for (int i = 0; i < 2; ++i)
  // {
  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  readInput(NULL);
  cmo_init_misc();

  CMO_PROF_FUNC_START();

  spectrum.updateJ = TRUE;

  getCPU(1, TIME_START, NULL);
  MULTIatmos(&atmos, &geometry);
  if (atmos.Stokes)
    Bproject();

  readAtomicModels();
  readMolecularModels();
  SortLambda();

  getBoundary(&geometry);

  // Background(write_analyze_output = TRUE, equilibria_only = FALSE);
  cmo_init_background();
  cmo_background(write_analyze_output = TRUE, equilibria_only = FALSE);
  convertScales(&atmos, &geometry);

  getProfiles();
  initSolution();
  initScatter();

  getCPU(1, TIME_POLL, "Total Initialize");

  /* --- Solve radiative transfer for active ingredients -- --------- */

  Iterate(input.NmaxIter, input.iterLimit);

  printf("Initial Stat Eq achieved. Waiting\n");
  getchar();

  double timeStep = 0.5;
  for (int i = 11;  i < 31;  i++) 
  {
    double di = (i - 20.0) / 3.0;
    atmos.T[i] *= 1.0 + 2.0 * exp(-SQ(di));
  }


  for (int t = 0; t < 3; ++t)
  {
    // Update the background
    cmo_background(FALSE, FALSE);
    // Really need to update collisions too -- this is in fact done by the background already. Woot woot.
    // Get the new Voigt profiles for the adjusted atmosphere
    getProfiles();

    printf("Time integration %d\n", t);
    time_dependent_iteration(timeStep, input.NmaxIter, input.iterLimit);
    printf("Iter %d done. Waiting\n", t);
    getchar();
  }




  CMO_PROF_FUNC_END();
  cmo_free_misc();
  // }
  // CMO: Removed temporarily
  // adjustStokesMode();
  // niter = 0;
  // while (niter < input.NmaxScatter) {
  //   if (solveSpectrum(FALSE, FALSE) <= input.iterLimit)
  //     break;
  //   niter++;
  // }
  /* --- Write output files --                         -------------- */

  // if (atmos.hydrostatic) {
  //   geometry.scale = COLUMN_MASS;
  //   convertScales(&atmos, &geometry);
  // }
  // getCPU(1, TIME_START, NULL);

  // writeInput();
  // writeAtmos(&atmos);
  // writeGeometry(&geometry);
  // writeSpectrum(&spectrum);
  // writeFlux(FLUX_DOT_OUT);

  // for (nact = 0; nact < atmos.Nactiveatom; nact++) {
  //   atom = atmos.activeatoms[nact];

  //   writeAtom(atom);
  //   writePopulations(atom);
  //   writeRadRate(atom);
  //   writeCollisionRate(atom);
  //   writeDamping(atom);
  // }
  // for (nact = 0; nact < atmos.Nactivemol; nact++) {
  //   molecule = atmos.activemols[nact];
  //   writeMolPops(molecule);
  // }

  // writeOpacity();

  // getCPU(1, TIME_POLL, "Write output");
}
#endif
/* ------- end ---------------------------- rhf1d.c ----------------- */
