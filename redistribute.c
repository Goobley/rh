/* ------- file: -------------------------- redistribute.c ----------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Apr  1 14:01:22 2009 --

       --------------------------                      ----------RH-- */

/* --- Administers iterations to redistribute intensity in PRD line while
       keeping the population number fixed. --         -------------- */

#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "accelerate.h"
#include "error.h"
#include "inputs.h"
#include "CmoProfile.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[];
static void scatter_sched(void* userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId);

/* ------- begin -------------------------- Redistribute.c ---------- */

void Redistribute(int NmaxIter, double iterLimit) {
  const char routineName[] = "Redistribute";
  register int kr, nact;

  bool_t quiet, accel, eval_operator, redistribute;
  enum Interpolation representation;
  int niter, Nlamu;
  double drho, drhomax, drhomaxa;
  Atom *atom;
  AtomicLine *line;
  CMO_PROF_FUNC_START();

  int nPrdLines = 0;
  for (nact = 0; nact < atmos.Nactiveatom; ++nact)
  {
    atom = atmos.activeatoms[nact];
    nPrdLines += atom->Nprd;
  }

  AtomicLine* prdLines[nPrdLines];
  int nPrd = 0;

  CMO_PROF_REGION_START("Reg: NgInit");
  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];

    /* --- Initialize structures for Ng acceleration PRD iteration -- */

    for (kr = 0; kr < atom->Nline; kr++) {
      line = &atom->line[kr];
      if (line->PRD)
      {
        prdLines[nPrd++] = line;
      }
      if (line->PRD && line->Ng_prd == NULL) {
        if (input.PRD_angle_dep == PRD_ANGLE_DEP)
          Nlamu = 2 * atmos.Nrays * line->Nlambda * atmos.Nspace;
        else
          Nlamu = line->Nlambda * atmos.Nspace;

        line->Ng_prd = NgInit(Nlamu, input.PRD_Ngdelay, input.PRD_Ngorder,
                              input.PRD_Ngperiod, line->rho_prd[0]);
      }
    }
  }
  CMO_PROF_REGION_END("Reg: NgInit");
  /* --- Iterate over scattering integral while keeping populations
         fixed --                                      -------------- */

  CMO_PROF_REGION_START("Reg: PRDIterate");
  niter = 1;
  // NOTE(cmo): Don't change the representation away from LINEAR, or non thread-safe itnerpolations are used.
  while (niter <= NmaxIter) {

    drhomaxa = 0.0;
    if (FALSE)
    {
    for (nact = 0; nact < atmos.Nactiveatom; nact++) {
      atom = atmos.activeatoms[nact];

      drhomax = 0.0;
      for (kr = 0; kr < atom->Nline; kr++) {
        line = &atom->line[kr];
        if (line->PRD) {
          switch (input.PRD_angle_dep) {
          case PRD_ANGLE_INDEP:
            PRDScatter(line, representation = LINEAR);
            break;

          case PRD_ANGLE_APPROX:
            PRDAngleApproxScatter(line, representation = LINEAR);
            break;

          case PRD_ANGLE_DEP:
            PRDAngleScatter(line, representation = LINEAR);
            break;
          }

          accel = Accelerate(line->Ng_prd, line->rho_prd[0]);
          sprintf(messageStr, "  PRD: iter #%d, atom %s, line %d,",
                  line->Ng_prd->count - 1, atom->ID, kr);
          drho = MaxChange(line->Ng_prd, messageStr, quiet = FALSE);
          sprintf(messageStr, (accel) ? " (accelerated)\n" : "\n");
          Error(MESSAGE, routineName, messageStr);

          drhomax = MAX(drho, drhomax);
        }
        drhomaxa = MAX(drhomax, drhomaxa);
      }
    }
    }
    else
    {
      CMO_PROF_REGION_START("REG: SCATTER");
      {
        struct sched_task task;
        scheduler_add(&input.sched, &task, scatter_sched, prdLines, nPrd, 1);
        scheduler_join(&input.sched, &task);
      }
      CMO_PROF_REGION_END("REG: SCATTER");
      CMO_PROF_REGION_START("REG: ACCEL");
      for (int krp = 0; krp < nPrd; ++krp)
      {
          line = prdLines[krp];
          atom = line->atom;
          kr = line->nLine;
          accel = Accelerate(line->Ng_prd, line->rho_prd[0]);
          sprintf(messageStr, "  PRD: iter #%d, atom %s, line %d,",
                  line->Ng_prd->count - 1, atom->ID, kr);
          drho = MaxChange(line->Ng_prd, messageStr, quiet = FALSE);
          sprintf(messageStr, (accel) ? " (accelerated)\n" : "\n");
          Error(MESSAGE, routineName, messageStr);

          drhomaxa = MAX(drho, drhomaxa);

      }
      CMO_PROF_REGION_END("REG: ACCEL");
    }
    
    /* --- Solve transfer equation with fixed populations -- -------- */

    // solveSpectrum(eval_operator = FALSE, redistribute = TRUE);
    solve_spectrum_redist(eval_operator = FALSE);

    if (drhomaxa < iterLimit)
      break;
    niter++;
  }
  CMO_PROF_REGION_END("Reg: PRDIterate");
  CMO_PROF_FUNC_END();
}
/* ------- end ---------------------------- Redistribute.c ---------- */

static void scatter_sched(void* userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId)
{
  // NOTE(cmo): Don't change the representation away from LINEAR, or non thread-safe itnerpolations are used.
  enum Interpolation representation;
  CMO_PROF_FUNC_START();
  AtomicLine** lines = (AtomicLine**)userdata;
  for (int l = range.start; l < range.end; ++l)
  {
    AtomicLine* line = lines[l];
    switch (input.PRD_angle_dep) {
    case PRD_ANGLE_INDEP:
      PRDScatter(line, representation = LINEAR);
      break;

    case PRD_ANGLE_APPROX:
      PRDAngleApproxScatter(line, representation = LINEAR);
      break;

    case PRD_ANGLE_DEP:
      PRDAngleScatter(line, representation = LINEAR);
      break;
    }

  }
  CMO_PROF_FUNC_END();
}

// static void accelerate_sched(void* userdata, struct scheduler *s, struct sched_task_partition range, sched_uint threadId)
// {
//   enum Interpolation representation;
//   CMO_PROF_FUNC_START();
//   AtomicLine* lines = (AtomicLine*)userdata;
//   for (int l = range.start; l < range.end; ++l)
//   {
//     AtomicLine* line = &lines[l];
//     accel = Accelerate(line->Ng_prd, line->rho_prd[0]);
//   }
//   CMO_PROF_FUNC_END();
// }
