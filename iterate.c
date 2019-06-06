/* ------- file: -------------------------- iterate.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek  (huitenbroek@nso.edu)
       Last modified: Tue Nov 16 15:31:48 2010 --

       --------------------------                      ----------RH-- */

/* --- Main iteration routine --                       -------------- */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "accelerate.h"
#include "error.h"
#include "statistics.h"
#include "inputs.h"
#include "CmoProfile.h"

typedef struct {
  bool_t eval_operator, redistribute;
  int nspect;
  double dJ;
} threadinfo;

typedef struct ThreadInfoRange {
  bool_t eval_operator, redistribute;
  int nspectStart, nspectEnd, lambdaMax;
  double dJ;
} ThreadInfoRange;

typedef struct FormalThreadInfo
{
  bool_t eval_operator, redistribute;
  int lambdaMax;
  double dJ;
} FormalThreadInfo;

typedef struct FormalThreadInfoPrd
{
  bool_t eval_operator, redistribute;
  int lambdaMax;
  double dJ;
  int* prdSpectralPoints;
} FormalThreadInfoPrd;

/* --- Function prototypes --                          -------------- */

void *Formal_pthread(void *argument);
void *Formal_pthread_range(void *argument);
static void Formal_sched(void *userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId);
static void Formal_complete_sched(void *userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId);
static void Formal_redist_sched(void *userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId);
static void accumulate_Gamma_sched(void *userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId);
static void accumulate_rates_lines_sched(void *userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId);
static void accumulate_rates_cont_sched(void *userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId);

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];

/* ------- begin -------------------------- Iterate.c --------------- */

void Iterate(int NmaxIter, double iterLimit) {
  const char routineName[] = "Iterate";
  register int niter, nact;
  double cswitch;

  bool_t eval_operator, write_analyze_output, equilibria_only;
  int Ngorder;
  double dpopsmax, PRDiterlimit;
  Atom *atom;
  Molecule *molecule;
  AtomicLine *line;        // Tiago: DELETE
  int i, mu, to_obs, lamu; // Tiago: DELETE
  CMO_PROF_FUNC_START();

  if (NmaxIter <= 0)
    return;
  getCPU(1, TIME_START, NULL);


  /* --- Initialize structures for Ng acceleration of population
         convergence --                                  ------------ */

  CMO_PROF_REGION_START("Reg: NgInit");
  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];
    atom->Ng_n = NgInit(atom->Nlevel * atmos.Nspace, input.Ngdelay,
                        input.Ngorder, input.Ngperiod, atom->n[0]);
  }
  for (nact = 0; nact < atmos.Nactivemol; nact++) {
    molecule = atmos.activemols[nact];
    Ngorder = (input.accelerate_mols) ? input.Ngorder : 0;

    molecule->Ng_nv = NgInit(molecule->Nv * atmos.Nspace, input.Ngdelay,
                             Ngorder, input.Ngperiod, molecule->nv[0]);
  }
  CMO_PROF_REGION_END("Reg: NgInit");
  /* --- Start of the main iteration loop --             ------------ */

  niter = 1;

  /* Collisional-radiative switching ? */
  if (input.crsw != 0.0)
    cswitch = input.crsw_ini;
  else
    cswitch = 1.0;

  /* PRD switching ? */
  if (input.prdsw > 0.0)
    input.prdswitch = 0.0;
  else
    input.prdswitch = 1.0;

  // while (niter <= NmaxIter && !StopRequested()) {
  while (niter <= NmaxIter) {
    CMO_PROF_REGION_START("Reg: Iterate");
    getCPU(2, TIME_START, NULL);

    CMO_PROF_REGION_START("Reg: Iterate: initGamma");
    for (nact = 0; nact < atmos.Nactiveatom; nact++)
      initGammaAtom(atmos.activeatoms[nact], cswitch);
    for (nact = 0; nact < atmos.Nactivemol; nact++)
      initGammaMolecule(atmos.activemols[nact]);
    CMO_PROF_REGION_END("Reg: Iterate: initGamma");

    /* --- Formal solution for all wavelengths --      -------------- */

    // solveSpectrum(eval_operator = TRUE, FALSE);
    solve_spectrum_complete(eval_operator = TRUE);

    /* --- Solve statistical equilibrium equations --  -------------- */

    sprintf(messageStr,
            "\n -- Iteration %3d, switch = %.2f, prd switch = %.2f\n", niter,
            cswitch, input.prdswitch);
    Error(MESSAGE, routineName, messageStr);
    dpopsmax = updatePopulations(niter);

    if (atmos.NPRDactive > 0) {

      /* --- Redistribute intensity in PRD lines if necessary -- ---- */

      if (input.PRDiterLimit < 0.0)
        PRDiterlimit = MAX(dpopsmax, -input.PRDiterLimit);
      else
        PRDiterlimit = input.PRDiterLimit;
      Redistribute(input.PRD_NmaxIter, PRDiterlimit);
    }

    sprintf(messageStr, "Total Iteration %3d", niter);
    getCPU(2, TIME_POLL, messageStr);

    if (dpopsmax < iterLimit && cswitch <= 1.0 && input.prdswitch == 1.0)
      break;
    niter++;

    if (input.solve_ne == ITERATION)
      // Background(write_analyze_output = TRUE, equilibria_only = FALSE);
      cmo_background(write_analyze_output = TRUE, equilibria_only = FALSE);

    /* Update collisional multiplier factor */
    if (input.crsw > 0)
      cswitch = MAX(1.0, cswitch * pow(0.1, 1. / input.crsw));

    /* Update PRD switching */
    if (input.prdsw > 0.0)
      input.prdswitch =
          MIN(1.0, input.prdsw * (double)(niter * niter)); // quadratic, for now

    if (atmos.hydrostatic) {
      if (!atmos.atoms[0].active) {
        sprintf(messageStr, "Can only perform hydrostatic equilibrium"
                            " for hydrogen active");
        Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      Hydrostatic(N_MAX_HSE_ITER, HSE_ITER_LIMIT);
    }
    CMO_PROF_REGION_END("Reg: Iterate");
  }

  // Tiago: temporary printouts to get PRD rho after iteration
  /*
   atom = atmos.activeatoms[0];
   line = &atom->line[0];

   switch (input.PRD_angle_dep) {
    case PRD_ANGLE_INDEP:
      printf("rho_prd = \n");
      for (i = 0; i < line->Nlambda; i++) {
        printf("%8.4f   %e   %e   %e   %e   %e\n", line->lambda[i],
   line->rho_prd[i][105], line->rho_prd[i][110], line->rho_prd[i][120],
   line->rho_prd[i][150], line->rho_prd[i][155]);
      }
      //exit(1);
      break;

    case PRD_ANGLE_DEP:
        for (mu = 0; mu < atmos.Nrays; mu++) {
          for (to_obs = 0; to_obs <= 1; to_obs++) {
           for (i = 0; i < line->Nlambda; i++) {
            lamu = 2*(atmos.Nrays*i + mu) + to_obs;
            if ((to_obs == 1) && (mu == 4))
            printf("%8.4f  %e   %e   %e   %e   %e\n", line->lambda[i],
   line->rho_prd[lamu][105],line->rho_prd[lamu][110],line->rho_prd[lamu][120],
   line->rho_prd[lamu][150], line->rho_prd[lamu][155] );
          }
        }
      }
      //exit(1);
      break;
   }
  */

  CMO_PROF_REGION_START("Reg: CleanUp");
  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];
    freeMatrix((void **)atom->Gamma);
    NgFree(atom->Ng_n);
  }
  for (nact = 0; nact < atmos.Nactivemol; nact++) {
    molecule = atmos.activemols[nact];
    freeMatrix((void **)molecule->Gamma);
    NgFree(molecule->Ng_nv);
  }
  CMO_PROF_REGION_END("Reg: CleanUp");

  getCPU(1, TIME_POLL, "Iteration Total");
  CMO_PROF_FUNC_END();
}
/* ------- end ---------------------------- Iterate.c --------------- */

/* ------- begin -------------------------- solveSpectrum.c --------- */

double solveSpectrum(bool_t eval_operator, bool_t redistribute) {
  register int nspect, n, nt, k;

  int Nthreads, lambda_max;
  double dJ, dJmax;
  pthread_t *thread_id;
  ThreadInfoRange *ti;
  /* --- Administers the formal solution for each wavelength. When
         input.Nthreads > 1 the solutions are performed concurrently
         in Nthreads threads. These are POSIX style threads.

    See: - David R. Butenhof, Programming with POSIX threads,
           Addison & Wesley.

         - Multithreaded Programming Guide, http://sun.docs.com
           (search for POSIX threads).

         When solveSpectrum is called with redistribute == TRUE only
         wavelengths that have an active PRD line are solved. The
         redistribute key is passed to the addtoRates routine via
         Formal so that only the radiative rates of PRD lines are
         updated. These are needed for the emission profile ratio \rho.
         --                                            -------------- */

  CMO_PROF_FUNC_START();
  getCPU(3, TIME_START, NULL);

  /* --- First zero the radiative rates --             -------------- */

  zeroRates(redistribute);
  lambda_max = 0;
  dJmax = 0.0;

  /* zero out J in gas parcel's frame */
  if (spectrum.updateJ && input.PRD_angle_dep == PRD_ANGLE_APPROX &&
      atmos.Nrays > 1 && atmos.NPRDactive > 0) {
    for (k = 0; k < atmos.Nspace; k++) {
      for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
        spectrum.Jgas[nspect][k] = 0.0;
      }
    }
  }

  if (input.Nthreads > 1) {

    if (FALSE)
    {

    /* --- If input.Nthreads positive then solve Nthreads wavelengths
           concurrently in separate threads --         -------------- */

    ti = (ThreadInfoRange *)malloc(input.Nthreads * sizeof(ThreadInfoRange));
    for (nt = 0; nt < input.Nthreads; nt++) {
      ti[nt].eval_operator = eval_operator;
      ti[nt].redistribute = redistribute;
    }
    thread_id = (pthread_t *)malloc(input.Nthreads * sizeof(pthread_t));

    /* --- Thread management is very simple. Submit a batch of as many
           as input.Nthreads at the same time, then wait till all of
           these have finished. There is no check on successful
           submission nor completion. --               -------------- */

    for (nspect = 0; nspect < spectrum.Nspect; nspect += input.Nthreads * input.workSize) {

      if (nspect + input.Nthreads * input.workSize <= spectrum.Nspect) {
        Nthreads = input.Nthreads;
      } else {
        Nthreads = (spectrum.Nspect - nspect) / input.workSize;
        if (((spectrum.Nspect - nspect) % input.workSize) != 0) {
          // Handle remaining odd-sized batch
          Nthreads += 1;
        }
      }
      /* --- Start batch of concurrent threads --      -------------- */

      for (nt = 0; nt < Nthreads; nt++) {
        ti[nt].nspectStart = nspect + nt * input.workSize;
        ti[nt].nspectEnd = nspect + (nt + 1) * input.workSize;
      }
      if (ti[Nthreads-1].nspectEnd > spectrum.Nspect) {
        ti[Nthreads-1].nspectEnd = spectrum.Nspect;
      }
      for (nt = 0; nt < Nthreads; nt++) {
        pthread_create(&thread_id[nt], &input.thread_attr, Formal_pthread_range, &ti[nt]);
      }
      /* --- Let the finished threads of the batch join again -- ---- */

      for (nt = 0; nt < Nthreads; nt++) {
        if (thread_id[nt]) {
          pthread_join(thread_id[nt], NULL);
          if (ti[nt].dJ > dJmax) {
            dJmax = ti[nt].dJ;
            lambda_max = ti[nt].lambdaMax;
          }
        }
      }
    }
    free(thread_id);
    free(ti);
    for (int nact = 0; nact < atmos.Nactiveatom; ++nact)
    {
      accumulate_Gamma(atmos.activeatoms[nact]);
      accumulate_rates_lines(atmos.activeatoms[nact]);
      accumulate_rates_cont(atmos.activeatoms[nact]);
    }

    }
    else
    {

      FormalThreadInfo* threadData = (FormalThreadInfo *)malloc(input.Nthreads * sizeof(FormalThreadInfo));
      for (int n = 0; n < input.Nthreads; ++n)
      {
        threadData[n].eval_operator = eval_operator;
        threadData[n].redistribute = redistribute;
        threadData[n].lambdaMax = 0;
        threadData[n].dJ = 0;
      }
      {
        struct sched_task task;
        scheduler_add(&input.sched, &task, Formal_sched, threadData, spectrum.Nspect, input.workSize);
        scheduler_join(&input.sched, &task);
      }

      for (int n = 0; n < input.Nthreads; ++n)
      {
	// THIS IS PRESUMABLY NOT ACCUMULATED PROPERLY, WHY?
        if (threadData[n].dJ > dJmax)
        {
          dJmax = threadData[n].dJ;
          lambda_max = threadData[n].lambdaMax;
        }
      }
      free(threadData);

      {
        struct sched_task tGamma;
        struct sched_task tRatesLines;
        struct sched_task tRatesCont;
        scheduler_add(&input.sched, &tGamma, accumulate_Gamma_sched, NULL, atmos.Nactiveatom, 1);
        scheduler_add(&input.sched, &tRatesLines, accumulate_rates_lines_sched, NULL, atmos.Nactiveatom, 1);
        scheduler_add(&input.sched, &tRatesCont, accumulate_rates_cont_sched, NULL, atmos.Nactiveatom, 1);
        scheduler_join(&input.sched, &tGamma);
        scheduler_join(&input.sched, &tRatesLines);
        scheduler_join(&input.sched, &tRatesCont);
      }
    }
    
  } else {

    /* --- Else call the solution for wavelengths sequentially -- --- */

    for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
      if (!redistribute ||
          (redistribute && containsPRDline(&spectrum.as[nspect]))) {
        dJ = Formal(nspect, eval_operator, redistribute, 0);
        if (dJ > dJmax) {
          dJmax = dJ;
          lambda_max = nspect;
        }
      }
    }
    for (int nact = 0; nact < atmos.Nactiveatom; ++nact)
    {
      accumulate_Gamma(atmos.activeatoms[nact]);
      accumulate_rates_lines(atmos.activeatoms[nact]);
      accumulate_rates_cont(atmos.activeatoms[nact]);
    }
  }

  sprintf(messageStr, " Spectrum max delta J = %6.4E (lambda#: %d)\n", dJmax,
          lambda_max);
  Error(MESSAGE, NULL, messageStr);

  getCPU(3, TIME_POLL,
         (eval_operator) ? "Spectrum & Operator" : "Solve Spectrum");

  CMO_PROF_FUNC_END();
  return dJmax;
}



double solve_spectrum_complete(bool_t eval_operator) {
  register int nspect, n, nt, k;

  int Nthreads, lambda_max;
  double dJ, dJmax;
  pthread_t *thread_id;
  ThreadInfoRange *ti;
  bool_t redistribute = FALSE;
  /* --- Administers the formal solution for each wavelength. When
         input.Nthreads > 1 the solutions are performed concurrently
         in Nthreads threads. These are POSIX style threads.

    See: - David R. Butenhof, Programming with POSIX threads,
           Addison & Wesley.

         - Multithreaded Programming Guide, http://sun.docs.com
           (search for POSIX threads).

         When solveSpectrum is called with redistribute == TRUE only
         wavelengths that have an active PRD line are solved. The
         redistribute key is passed to the addtoRates routine via
         Formal so that only the radiative rates of PRD lines are
         updated. These are needed for the emission profile ratio \rho.
         --                                            -------------- */

  CMO_PROF_FUNC_START();
  getCPU(3, TIME_START, NULL);

  /* --- First zero the radiative rates --             -------------- */

  zeroRates(redistribute);
  lambda_max = 0;
  dJmax = 0.0;

  /* zero out J in gas parcel's frame */
  if (spectrum.updateJ && input.PRD_angle_dep == PRD_ANGLE_APPROX &&
      atmos.Nrays > 1 && atmos.NPRDactive > 0) {
    for (k = 0; k < atmos.Nspace; k++) {
      for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
        spectrum.Jgas[nspect][k] = 0.0;
      }
    }
  }

  if (input.Nthreads > 1) {
    FormalThreadInfo* threadData = (FormalThreadInfo *)malloc(input.Nthreads * sizeof(FormalThreadInfo));
    for (int n = 0; n < input.Nthreads; ++n)
    {
      threadData[n].eval_operator = eval_operator;
      threadData[n].redistribute = redistribute;
      threadData[n].lambdaMax = 0;
      threadData[n].dJ = 0;
    }
    {
      struct sched_task task;
      scheduler_add(&input.sched, &task, Formal_complete_sched, threadData, spectrum.Nspect, input.workSize);
      scheduler_join(&input.sched, &task);
    }

    for (int n = 0; n < input.Nthreads; ++n)
    {
      if (threadData[n].dJ > dJmax)
      {
        dJmax = threadData[n].dJ;
        lambda_max = threadData[n].lambdaMax;
      }
    }
    free(threadData);

    {
      struct sched_task tGamma;
      struct sched_task tRatesLines;
      struct sched_task tRatesCont;
      scheduler_add(&input.sched, &tGamma, accumulate_Gamma_sched, NULL, atmos.Nactiveatom, 1);
      scheduler_add(&input.sched, &tRatesLines, accumulate_rates_lines_sched, NULL, atmos.Nactiveatom, 1);
      scheduler_add(&input.sched, &tRatesCont, accumulate_rates_cont_sched, NULL, atmos.Nactiveatom, 1);
      scheduler_join(&input.sched, &tGamma);
      scheduler_join(&input.sched, &tRatesLines);
      scheduler_join(&input.sched, &tRatesCont);
    }
  
  } else {

    /* --- Else call the solution for wavelengths sequentially -- --- */

    for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
      if (!redistribute ||
          (redistribute && containsPRDline(&spectrum.as[nspect]))) {
        dJ = Formal(nspect, eval_operator, redistribute, 0);
        if (dJ > dJmax) {
          dJmax = dJ;
          lambda_max = nspect;
        }
      }
    }
    for (int nact = 0; nact < atmos.Nactiveatom; ++nact)
    {
      accumulate_Gamma(atmos.activeatoms[nact]);
      accumulate_rates_lines(atmos.activeatoms[nact]);
      accumulate_rates_cont(atmos.activeatoms[nact]);
    }
  }

  sprintf(messageStr, " Spectrum max delta J = %6.4E (lambda#: %d)\n", dJmax,
          lambda_max);
  Error(MESSAGE, NULL, messageStr);

  getCPU(3, TIME_POLL,
         (eval_operator) ? "Spectrum & Operator" : "Solve Spectrum");

  CMO_PROF_FUNC_END();
  return dJmax;
}



double solve_spectrum_redist(bool_t eval_operator) {
  register int nspect, n, nt, k;

  int Nthreads, lambda_max;
  double dJ, dJmax;
  pthread_t *thread_id;
  ThreadInfoRange *ti;

  bool_t redistribute = TRUE;
  /* --- Administers the formal solution for each wavelength. When
         input.Nthreads > 1 the solutions are performed concurrently
         in Nthreads threads. These are POSIX style threads.

    See: - David R. Butenhof, Programming with POSIX threads,
           Addison & Wesley.

         - Multithreaded Programming Guide, http://sun.docs.com
           (search for POSIX threads).

         When solveSpectrum is called with redistribute == TRUE only
         wavelengths that have an active PRD line are solved. The
         redistribute key is passed to the addtoRates routine via
         Formal so that only the radiative rates of PRD lines are
         updated. These are needed for the emission profile ratio \rho.
         --                                            -------------- */

  CMO_PROF_FUNC_START();
  getCPU(3, TIME_START, NULL);

  /* --- First zero the radiative rates --             -------------- */

  zeroRates(redistribute);
  lambda_max = 0;
  dJmax = 0.0;

  /* zero out J in gas parcel's frame */
  if (spectrum.updateJ && input.PRD_angle_dep == PRD_ANGLE_APPROX &&
      atmos.Nrays > 1 && atmos.NPRDactive > 0) {
    for (k = 0; k < atmos.Nspace; k++) {
      for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
        spectrum.Jgas[nspect][k] = 0.0;
      }
    }
  }

  if (input.Nthreads > 1) {

    int prdPoints[spectrum.Nspect];
    int prdIdx = 0;
    for (nspect = 0; nspect < spectrum.Nspect; ++nspect)
    {
      if (containsPRDline(&spectrum.as[nspect]))
      {
        prdPoints[prdIdx++] = nspect;
      }
    }

    // FormalThreadInfo* threadData = (FormalThreadInfo *)malloc(input.Nthreads * sizeof(FormalThreadInfo));
    FormalThreadInfoPrd threadData[input.Nthreads];
    for (int n = 0; n < input.Nthreads; ++n)
    {
      threadData[n].eval_operator = eval_operator;
      threadData[n].redistribute = redistribute;
      threadData[n].lambdaMax = 0;
      threadData[n].dJ = 0;
      threadData[n].prdSpectralPoints = prdPoints;
    }
    {
      struct sched_task task;
      scheduler_add(&input.sched, &task, Formal_redist_sched, threadData, prdIdx, input.workSize);
      scheduler_join(&input.sched, &task);
    }

    for (int n = 0; n < input.Nthreads; ++n)
    {
      if (threadData[n].dJ > dJmax)
      {
        dJmax = threadData[n].dJ;
        lambda_max = threadData[n].lambdaMax;
      }
    }
    // free(threadData);

    {
      // Nothing is ever added to the continuum rates in PRD
      struct sched_task tGamma;
      struct sched_task tRatesLines;
      // struct sched_task tRatesCont;
      scheduler_add(&input.sched, &tGamma, accumulate_Gamma_sched, NULL, atmos.Nactiveatom, 1);
      scheduler_add(&input.sched, &tRatesLines, accumulate_rates_lines_sched, NULL, atmos.Nactiveatom, 1);
      // scheduler_add(&input.sched, &tRatesCont, accumulate_rates_cont_sched, NULL, atmos.Nactiveatom, 1);
      scheduler_join(&input.sched, &tGamma);
      scheduler_join(&input.sched, &tRatesLines);
      // scheduler_join(&input.sched, &tRatesCont);
    }
    
  } else {

    /* --- Else call the solution for wavelengths sequentially -- --- */

    for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
      if (!redistribute ||
          (redistribute && containsPRDline(&spectrum.as[nspect]))) {
        dJ = Formal(nspect, eval_operator, redistribute, 0);
        if (dJ > dJmax) {
          dJmax = dJ;
          lambda_max = nspect;
        }
      }
    }
    for (int nact = 0; nact < atmos.Nactiveatom; ++nact)
    {
      accumulate_Gamma(atmos.activeatoms[nact]);
      accumulate_rates_lines(atmos.activeatoms[nact]);
      accumulate_rates_cont(atmos.activeatoms[nact]);
    }
  }

  sprintf(messageStr, " Spectrum max delta J = %6.4E (lambda#: %d)\n", dJmax,
          lambda_max);
  Error(MESSAGE, NULL, messageStr);

  getCPU(3, TIME_POLL,
         (eval_operator) ? "Spectrum & Operator" : "Solve Spectrum");

  CMO_PROF_FUNC_END();
  return dJmax;
}
/* ------- end ---------------------------- solveSpectrum.c --------- */

/* ------- begin -------------------------- Formal_pthread.c -------- */

void *Formal_pthread(void *argument) {
  threadinfo *ti = (threadinfo *)argument;
  int nt = (ti->nspect / input.workSize) % input.Nthreads;
  /* --- Threads wrapper around Formal --              -------------- */

  ti->dJ = Formal(ti->nspect, ti->eval_operator, ti->redistribute, nt);

  return (NULL);
}
/* ------- end ---------------------------- Formal_pthread.c -------- */

void *Formal_pthread_range(void *argument) {
  ThreadInfoRange *ti = (ThreadInfoRange *)argument;
  int nspect, lambdaMax;
  double dJ, dJMax;
  dJMax = 0.0;
  lambdaMax = 0;
  int nt = (ti->nspectStart / input.workSize) % input.Nthreads;

  /* --- Threads wrapper around Formal --              -------------- */
  for (nspect = ti->nspectStart; nspect < ti->nspectEnd; ++nspect) {
    if (!ti->redistribute ||
        (ti->redistribute && containsPRDline(&spectrum.as[nspect]))) {
      dJ = Formal(nspect, ti->eval_operator, ti->redistribute, nt);
      if (dJ > dJMax) {
        lambdaMax = nspect;
        dJMax = dJ;
      }
    }
  }

  ti->dJ = dJMax;
  ti->lambdaMax = lambdaMax;

  return (NULL);
}

static void Formal_sched(void *userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId)
{
  CMO_PROF_FUNC_START();
  FormalThreadInfo* ti = &((FormalThreadInfo *)userdata)[threadId];

  for (int nspect = range.start; nspect < range.end; ++nspect)
  {
    if (!ti->redistribute ||
        (ti->redistribute && containsPRDline(&spectrum.as[nspect])))
    {
      double dJ = Formal(nspect, ti->eval_operator, ti->redistribute, threadId);
      if (dJ > ti->dJ)
      {
        ti->dJ = dJ;
        ti->lambdaMax = nspect;
      }
    }
  }
  CMO_PROF_FUNC_END();
}

static void Formal_complete_sched(void *userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId)
{
  CMO_PROF_FUNC_START();
  FormalThreadInfo* ti = &((FormalThreadInfo *)userdata)[threadId];

  for (int nspect = range.start; nspect < range.end; ++nspect)
  {
    // if (!ti->redistribute ||
    //     (ti->redistribute && containsPRDline(&spectrum.as[nspect])))
    // {
      double dJ = Formal(nspect, ti->eval_operator, ti->redistribute, threadId);
      if (dJ > ti->dJ)
      {
        ti->dJ = dJ;
        ti->lambdaMax = nspect;
      // }
    }
  }
  CMO_PROF_FUNC_END();
}

static void Formal_redist_sched(void *userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId)
{
  CMO_PROF_FUNC_START();
  FormalThreadInfoPrd* ti = &((FormalThreadInfoPrd *)userdata)[threadId];

  for (int nspect = range.start; nspect < range.end; ++nspect)
  {
    // if (!ti->redistribute ||
    //     (ti->redistribute && containsPRDline(&spectrum.as[nspect])))
    // {
      double dJ = Formal(ti->prdSpectralPoints[nspect], ti->eval_operator, ti->redistribute, threadId);
      if (dJ > ti->dJ)
      {
        ti->dJ = dJ;
        ti->lambdaMax = ti->prdSpectralPoints[nspect];
      // }
    }
  }
  CMO_PROF_FUNC_END();
}

static void accumulate_Gamma_sched(void *userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId)
{
  CMO_PROF_FUNC_START();
  for (int nact = range.start; nact < range.end; ++nact)
  {
    accumulate_Gamma(atmos.activeatoms[nact]);
  }
  CMO_PROF_FUNC_END();
}

static void accumulate_rates_lines_sched(void *userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId)
{
  CMO_PROF_FUNC_START();
  for (int nact = range.start; nact < range.end; ++nact)
  {
    accumulate_rates_lines(atmos.activeatoms[nact]);
  }
  CMO_PROF_FUNC_END();
}

static void accumulate_rates_cont_sched(void *userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId)
{
  CMO_PROF_FUNC_START();
  for (int nact = range.start; nact < range.end; ++nact)
  {
    accumulate_rates_cont(atmos.activeatoms[nact]);
  }
  CMO_PROF_FUNC_END();
}

