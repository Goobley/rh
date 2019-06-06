/* ------- file: -------------------------- fillgamma.c -------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Jul 24 12:34:18 2009 --

       --------------------------                      ----------RH-- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <alloca.h>

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "atom.h"
#include "inputs.h"
#include "constant.h"

#include "CmoProfile.h"

/* --- Routines for wavelength- and angle-integrated contributions to
       the crosscoupling, Gamma matrix, and radiative rates.

       Convention: \Gamma_ij = Gamma[i][j] represents the
                   transition j --> i, so that \Gamma_ij * n_j
                   is the rate per sec out of level j to level i.
       --                                              -------------- */

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];

/* ------- begin -------------------------- initGammaAtom.c --------- */

void initGammaAtom(Atom *atom, double cswitch) {
  register int ij, k;

  /* --- Add the fixed rates into Gamma --             -------------- */

  for (ij = 0; ij < SQ(atom->Nlevel); ij++) {
    for (k = 0; k < atmos.Nspace; k++)
      atom->Gamma[ij][k] = atom->C[ij][k] * cswitch;
  }
}
/* ------- end ---------------------------- initGammaAtom.c --------- */

/* ------- begin -------------------------- initGammaMolecule.c ----- */

void initGammaMolecule(Molecule *molecule) {
  register int ij, ji, k, vi, vj;

  /* --- Add the fixed rates into Gamma --             -------------- */

  for (vi = 0; vi < molecule->Nv - 1; vi++) {
    vj = vi + 1;
    ij = vi * molecule->Nv + vj;
    ji = vj * molecule->Nv + vi;
    for (k = 0; k < atmos.Nspace; k++) {
      if (molecule->n[k]) {
        molecule->Gamma[ij][k] = molecule->C_ul[k];
        molecule->Gamma[ji][k] = molecule->C_ul[k] * molecule->nvstar[vj][k] /
                                 molecule->nvstar[vi][k];
      }
    }
  }
}
/* ------- end ---------------------------- initGammaMolecule.c ----- */

/* ------- begin -------------------------- addtoGamma.c ------------ */

void addtoGamma(int nspect, double wmu, double *I, double *Psi, int threadId) {
// void addtoGamma(int nspect, double wmu, double *I, double *Psi) {
  const char routineName[] = "addtoGamma";
  register int nact, n, k, m;

  int i, j, ij, ji, jp, nt;
  double twohnu3_c2, twohc, wlamu, *Ieff, *Stokes_Q, *Stokes_U, *Stokes_V,
      *eta_Q, *eta_U, *eta_V;

  Atom *atom;
  AtomicLine *line;
  AtomicContinuum *continuum;
  Molecule *molecule;
  MolecularLine *mrt;
  ActiveSet *as;

  CMO_PROF_FUNC_START();

  nt = threadId;

  twohc = 2.0 * HPLANCK * CLIGHT / CUBE(NM_TO_M);

  as = &spectrum.as[nspect];
  // // nt = nspect % input.Nthreads;
  // if (input.Nthreads > 1)
  // {
  //   nt = (nspect / input.workSize) % input.Nthreads;
  // } else {
  //   nt = nspect % input.Nthreads;
  // }

  if (containsActive(as)) {
    Ieff = (double *)malloc(atmos.Nspace * sizeof(double));

    if (input.StokesMode == FULL_STOKES && containsPolarized(as)) {

      /* --- Use pointers to the bottom 3/4 of I and
             atom->rhth.eta --                         -------------- */

      Stokes_Q = I + atmos.Nspace;
      Stokes_U = I + 2 * atmos.Nspace;
      Stokes_V = I + 3 * atmos.Nspace;
    }
  }
  /* --- Contributions from the active transitions in atoms -- ------ */

  double scr_ji[atmos.Nspace];
  double scr_ij[atmos.Nspace];
  double scr_crossCouple[atmos.Nspace];
  double scr_otherUpper[atmos.Nspace];

  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];

    if (as->Nactiveatomrt[nact] > 0) {
      if (input.StokesMode == FULL_STOKES && containsPolarized(as)) {

        eta_Q = atom->rhth[nt].eta + atmos.Nspace;
        eta_U = atom->rhth[nt].eta + 2 * atmos.Nspace;
        eta_V = atom->rhth[nt].eta + 3 * atmos.Nspace;

        for (k = 0; k < atmos.Nspace; k++) {
          Ieff[k] =
              I[k] + Stokes_Q[k] + Stokes_U[k] + Stokes_V[k] -
              Psi[k] * (atom->rhth[nt].eta[k] + eta_Q[k] + eta_U[k] + eta_V[k]);
        }
      } else {
        for (k = 0; k < atmos.Nspace; k++) {
          Ieff[k] = I[k] - Psi[k] * atom->rhth[nt].eta[k];
        }
      }
    }

    for (n = 0; n < as->Nactiveatomrt[nact]; n++) {
      switch (as->art[nact][n].type) {
      case ATOMIC_LINE:
        line = as->art[nact][n].ptype.line;
        i = line->i;
        j = line->j;
        twohnu3_c2 = line->Aji / line->Bji;
        break;

      case ATOMIC_CONTINUUM:
        continuum = as->art[nact][n].ptype.continuum;
        i = continuum->i;
        j = continuum->j;
        twohnu3_c2 = twohc / CUBE(spectrum.lambda[nspect]);
        break;

      default:
        sprintf(messageStr, "Invalid transition type");
        Error(ERROR_LEVEL_1, routineName, messageStr);
        twohnu3_c2 = 0.0;
      }

      ij = i * atom->Nlevel + j;
      ji = j * atom->Nlevel + i;
      for (k = 0; k < atmos.Nspace; k++) {
        wlamu = atom->rhth[nt].Vij[n][k] * atom->rhth[nt].wla[n][k] * wmu;
        scr_ji[k] = Ieff[k] * wlamu;
        scr_ij[k] = (twohnu3_c2 + Ieff[k]) * atom->rhth[nt].gij[n][k] * wlamu;
        scr_crossCouple[k] = atom->rhth[nt].chi_up[i][k] * Psi[k] * 
                             atom->rhth[nt].Uji_down[j][k]* wmu;
      }


      // if (input.Nthreads > 1)
      //   pthread_mutex_lock(&atom->Gamma_lock);


      for (k = 0; k < atmos.Nspace; k++) {
        atom->rhacc[nt].Gamma[ji][k] += scr_ji[k];
        atom->rhacc[nt].Gamma[ij][k] += scr_ij[k];
      }
      /* --- Cross-coupling terms, currently only for Stokes_I -- --- */

      for (k = 0; k < atmos.Nspace; k++) {
        atom->rhacc[nt].Gamma[ij][k] -= scr_crossCouple[k];
      }
      // if (input.Nthreads > 1)
      //   pthread_mutex_unlock(&atom->Gamma_lock);
      /* --- If rt->i is also an upper level of another transition that
             is active at this wavelength then Gamma[ji] needs to be
             updated as well --                        -------------- */

      // cmo: Look at optimising this inside the lock
      for (m = 0; m < as->Nactiveatomrt[nact]; m++) {
        switch (as->art[nact][m].type) {
        case ATOMIC_LINE:
          jp = as->art[nact][m].ptype.line->j;
          break;
        case ATOMIC_CONTINUUM:
          jp = as->art[nact][m].ptype.continuum->j;
          break;
        default:;
        }
        if (jp == i) {
          for (k = 0; k < atmos.Nspace; k++) {
            scr_otherUpper[k] = atom->rhth[nt].chi_down[j][k] * Psi[k] *
                                atom->rhth[nt].Uji_down[i][k] * wmu;
          }
          // if (input.Nthreads > 1)
          //   pthread_mutex_lock(&atom->Gamma_lock);
          for (k = 0; k < atmos.Nspace; k++) {
            atom->rhacc[nt].Gamma[ji][k] += scr_otherUpper[k];
          }
          // if (input.Nthreads > 1)
          //   pthread_mutex_unlock(&atom->Gamma_lock);
        }
      }
      // if (input.Nthreads > 1)
      //   pthread_mutex_unlock(&atom->Gamma_lock);
    }
  }
  /* --- Add the active molecular contributions --     -------------- */

  for (nact = 0; nact < atmos.Nactivemol; nact++) {
    molecule = atmos.activemols[nact];

    for (n = 0; n < as->Nactivemolrt[nact]; n++) {
      switch (as->mrt[nact][n].type) {

      case VIBRATION_ROTATION:
        mrt = as->mrt[nact][n].ptype.vrline;
        i = mrt->vi;
        j = mrt->vj;
        twohnu3_c2 = mrt->Aji / mrt->Bji;
        break;

      default:
        sprintf(messageStr, "Invalid transition type");
        Error(ERROR_LEVEL_1, routineName, messageStr);
        twohnu3_c2 = 0.0;
      }

      if (input.Nthreads > 1)
        pthread_mutex_lock(&molecule->Gamma_lock);

      /* --- In case of molecular vibration-rotation transitions -- - */

      ij = i * molecule->Nv + j;
      ji = j * molecule->Nv + i;

      for (k = 0; k < atmos.Nspace; k++) {
        if (molecule->n[k]) {
          wlamu =
              molecule->rhth[nt].Vij[n][k] * molecule->rhth[nt].wla[n][k] * wmu;
          molecule->Gamma[ji][k] += I[k] * wlamu;
          molecule->Gamma[ij][k] +=
              molecule->rhth[nt].gij[n][k] * (twohnu3_c2 + I[k]) * wlamu;
        }
      }
      if (input.Nthreads > 1)
        pthread_mutex_unlock(&molecule->Gamma_lock);
    }
  }

  if (containsActive(as))
    free(Ieff);
  CMO_PROF_FUNC_END();
}
/* ------- end ---------------------------- addtoGamma.c ------------ */

/* ------- begin -------------------------- addtoCoupling.c --------- */

void addtoCoupling(int nspect, int threadId) {
  const char routineName[] = "addtoCoupling";
  register int nact, n, k;

  int i, j, nt;
  double twohnu3_c2, chicc, twohc, *n_i, *n_j;
  Atom *atom;
  AtomicLine *line;
  AtomicContinuum *continuum;
  ActiveSet *as;

  CMO_PROF_FUNC_START();

  twohc = 2.0 * HPLANCK * CLIGHT / CUBE(NM_TO_M);

  as = &spectrum.as[nspect];
  // nt = nspect % input.Nthreads;
  nt = threadId;
  // if (input.Nthreads > 1)
  // {
  //   nt = (nspect / input.workSize) % input.Nthreads;
  // } else {
  //   nt = nspect % input.Nthreads;
  // }

  /* --- Zero the cross coupling matrices --           -------------- */

  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];

    for (n = 0; n < as->Nlower[nact]; n++) {
      i = as->lower_levels[nact][n];
      for (k = 0; k < atmos.Nspace; k++)
        atom->rhth[nt].chi_up[i][k] = 0.0;
    }
    for (n = 0; n < as->Nupper[nact]; n++) {
      j = as->upper_levels[nact][n];
      for (k = 0; k < atmos.Nspace; k++) {
        atom->rhth[nt].chi_down[j][k] = 0.0;
        atom->rhth[nt].Uji_down[j][k] = 0.0;
      }
    }
  }
  /* --- Gather terms for cross-coupling between overlapping
     transitions --                                -------------- */

  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];

    for (n = 0; n < as->Nactiveatomrt[nact]; n++) {
      switch (as->art[nact][n].type) {
      case ATOMIC_LINE:
        line = as->art[nact][n].ptype.line;
        i = line->i;
        j = line->j;
        n_i = atom->n[i];
        n_j = atom->n[j];
        twohnu3_c2 = line->Aji / line->Bji;
        break;

      case ATOMIC_CONTINUUM:
        continuum = as->art[nact][n].ptype.continuum;
        i = continuum->i;
        j = continuum->j;
        n_i = atom->n[i];
        n_j = atom->n[j];
        twohnu3_c2 = twohc / CUBE(spectrum.lambda[nspect]);
        break;

      default:
        sprintf(messageStr, "Invalid transition type");
        Error(ERROR_LEVEL_1, routineName, messageStr);
        twohnu3_c2 = 0.0;
      }
      /* --- Evaluate the cross-coupling coefficients -- ------------ */

      if (twohnu3_c2) {
        for (k = 0; k < atmos.Nspace; k++) {
          chicc = atom->rhth[nt].Vij[n][k] * atom->rhth[nt].wla[n][k] *
                  (n_i[k] - atom->rhth[nt].gij[n][k] * n_j[k]);
          atom->rhth[nt].chi_up[i][k] += chicc;
          atom->rhth[nt].chi_down[j][k] += chicc;

          atom->rhth[nt].Uji_down[j][k] +=
              twohnu3_c2 * atom->rhth[nt].gij[n][k] * atom->rhth[nt].Vij[n][k];
        }
      }
    }
  }
  CMO_PROF_FUNC_END();
}
/* ------- end ---------------------------- addtoCoupling.c --------- */

/* ------- begin -------------------------- zeroRates.c ------------- */

void zeroRates(bool_t redistribute) {
  register int kr, k, n;

  Atom *atom;

  /* --- Initialize the radiative rates for atomic transitions.

         When redistribute == TRUE only the rates of PRD lines are
         initialized.
         --                                            -------------- */

  CMO_PROF_FUNC_START();
  for (n = 0; n < atmos.Natom; n++) {
    atom = &atmos.atoms[n];
    if (atom->active) {
      for (kr = 0; kr < atom->Nline; kr++) {
        if (!redistribute || (redistribute && atom->line[kr].PRD)) {
          for (k = 0; k < atmos.Nspace; k++) {
            atom->line[kr].Rij[k] = 0.0;
            atom->line[kr].Rji[k] = 0.0;
          }
        }
      }
      if (!redistribute) {
        for (kr = 0; kr < atom->Ncont; kr++) {
          for (k = 0; k < atmos.Nspace; k++) {
            atom->continuum[kr].Rij[k] = 0.0;
            atom->continuum[kr].Rji[k] = 0.0;
          }
        }
      }
    }
  }
  CMO_PROF_FUNC_END();
}
/* ------- end ---------------------------- zeroRates.c ------------- */

/* ------- begin -------------------------- addtoRates.c ------------ */

void addtoRates(int nspect, int mu, bool_t to_obs, double wmu, double *I,
                bool_t redistribute, int threadId) {
  register int nact, n, k;

  int la, lamu, nt, nTrans;
  double twohnu3_c2, twohc, hc_4PI, Bijxhc_4PI, wlamu, *Rij, *Rji, up_rate,
      *Stokes_Q, *Stokes_U, *Stokes_V;

  ActiveSet *as;
  Atom *atom;
  AtomicLine *line;
  AtomicContinuum *continuum;
  pthread_mutex_t *rate_lock;
  
  CMO_PROF_FUNC_START();

  /* --- Calculate the radiative rates for atomic transitions.

         When redistribute == TRUE only radiative rates of PRD lines
         are evaluated. These are needed in the iterative update of the
         emission profile ratio \rho.
         --                                            -------------- */

  twohc = 2.0 * HPLANCK * CLIGHT / CUBE(NM_TO_M);

  as = &spectrum.as[nspect];
  nt = threadId;
  // if (input.Nthreads > 1)
  // {
  //   nt = (nspect / input.workSize) % input.Nthreads;
  // } else {
  //   nt = nspect % input.Nthreads;
  // }


  if (input.StokesMode == FULL_STOKES && containsPolarized(as)) {

    /* --- Use pointers to the bottom 3/4 of I and as->eta -- ------- */

    Stokes_Q = I + atmos.Nspace;
    Stokes_U = I + 2 * atmos.Nspace;
    Stokes_V = I + 3 * atmos.Nspace;
  }

  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];

    for (n = 0; n < as->Nactiveatomrt[nact]; n++) {
      switch (as->art[nact][n].type) {
      case ATOMIC_LINE:
        line = as->art[nact][n].ptype.line;
        if (redistribute && !line->PRD)
          Rij = NULL;
        else {
          nTrans = line->nLine;
          atom->rhacc[nt].lineRatesDirty[nTrans] = TRUE;
          Rij = atom->rhacc[nt].RijLine[nTrans];
          Rji = atom->rhacc[nt].RjiLine[nTrans];
          twohnu3_c2 = line->Aji / line->Bji;

          rate_lock = &line->rate_lock;
        }
        break;

      case ATOMIC_CONTINUUM:
        if (redistribute)
          Rij = NULL;
        else {
          continuum = as->art[nact][n].ptype.continuum;
          nTrans = continuum->nCont;
          Rij = atom->rhacc[nt].RijCont[nTrans];
          Rji = atom->rhacc[nt].RjiCont[nTrans];
          twohnu3_c2 = twohc / CUBE(spectrum.lambda[nspect]);

          rate_lock = &continuum->rate_lock;
        }
        break;

      default:
        Rij = NULL;
      }
      /* --- Convention: Rij is the rate for transition i -> j -- ----- */

      if (Rij != NULL) {
        double scr_Rij[atmos.Nspace];
        double scr_Rji[atmos.Nspace];
        for (k = 0; k < atmos.Nspace; k++) {
          wlamu = atom->rhth[nt].Vij[n][k] * atom->rhth[nt].wla[n][k] * wmu;
          scr_Rij[k] = I[k] * wlamu;
          scr_Rji[k] = atom->rhth[nt].gij[n][k] * (twohnu3_c2 + I[k]) * wlamu;
        }
        // if (input.Nthreads > 1)
        //   pthread_mutex_lock(rate_lock);

        for (k = 0; k < atmos.Nspace; k++) {
          Rij[k] += scr_Rij[k];
          Rji[k] += scr_Rji[k];
        }

        // if (input.Nthreads > 1)
        //   pthread_mutex_unlock(rate_lock);
      }
    }
  }
  CMO_PROF_FUNC_END();
}
/* ------- end ---------------------------- addtoRates.c ------------ */

void accumulate_Gamma(Atom* atom)
{
  // printf("Called %s\n", (atom->rhacc[7].Gamma == NULL) ? "NULL" : "FINE");
  // printf("%s\n", atom->ID);
  CMO_PROF_FUNC_START();
  for (int ij = 0; ij < SQ(atom->Nlevel); ++ij)
  {
    for (int k = 0; k < atmos.Nspace; ++k)
    {
      for (int nt = 0; nt < input.Nthreads; ++nt)
      {
        atom->Gamma[ij][k] += atom->rhacc[nt].Gamma[ij][k];
        atom->rhacc[nt].Gamma[ij][k] = 0.0;
      }
    }
  }
  CMO_PROF_FUNC_END();
}

void accumulate_rates_lines(Atom* atom)
{
  CMO_PROF_FUNC_START();
  for (int kr = 0; kr < atom->Nline; ++kr)
  {
    bool_t dirty = FALSE;
    for (int nt = 0; nt < input.Nthreads; ++nt)
    {
      if (atom->rhacc[nt].lineRatesDirty[kr])
      {
        dirty = TRUE;
      }
      atom->rhacc[nt].lineRatesDirty[kr] = FALSE;
    }
    if (dirty)
    {
      for (int k = 0; k < atmos.Nspace; ++k)
      {
        for (int nt = 0; nt < input.Nthreads; ++nt)
        {
          atom->line[kr].Rij[k] += atom->rhacc[nt].RijLine[kr][k];
          atom->line[kr].Rji[k] += atom->rhacc[nt].RjiLine[kr][k];
          atom->rhacc[nt].RijLine[kr][k] = 0.0;
          atom->rhacc[nt].RjiLine[kr][k] = 0.0;
        }
      }
    }
  }
  CMO_PROF_FUNC_END();
}

void accumulate_rates_cont(Atom* atom)
{
  CMO_PROF_FUNC_START();
  for (int kr = 0; kr < atom->Ncont; ++kr)
  {
    for (int k = 0; k < atmos.Nspace; ++k)
    {
      for (int nt = 0; nt < input.Nthreads; ++nt)
      {
        atom->continuum[kr].Rij[k] += atom->rhacc[nt].RijCont[kr][k];
        atom->continuum[kr].Rji[k] += atom->rhacc[nt].RjiCont[kr][k];
        atom->rhacc[nt].RijCont[kr][k] = 0.0;
        atom->rhacc[nt].RjiCont[kr][k] = 0.0;
      }
    }
  }
  CMO_PROF_FUNC_END();
}
