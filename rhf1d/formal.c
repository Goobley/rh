/* ------- file: -------------------------- formal.c ----------------

       Version:       rh1.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Jan 14 16:14:23 2010 --

       --------------------------                      ----------RH-- */

/* --- Formal solution with given source function, and allowing for
       a PRD emission profile, polarized line radiation, polarized
       background lines, and scattering background polarization -- -- */

#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "constant.h"
#include "background.h"
#include "inputs.h"
#include "error.h"
#include "xdr.h"
// #define CMO_NO_PROF
#include "CmoProfile.h"

/* --- Function prototypes --                          -------------- */

void loadBackground(int la, int mu, bool_t to_obs);
void cmo_load_background(int la, int mu, bool_t to_obs);

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];

double Formal_sc_no_pol_cmo(int nspect, bool_t eval_operator, bool_t redistribute, int threadId) {
  const char routineName[] = "Formal";
  register int k, mu, n;

  bool_t initialize, boundbound, polarized_as, polarized_c, PRD_angle_dep,
      to_obs, solveStokes, angle_dep;
  enum FeautrierOrder F_order;
  int Nrays = atmos.Nrays, lamuk;
  long Nspace = atmos.Nspace;
  long int idx, idx0;
  double *I, *chi, *S, **Ipol, **Spol, *Psi, *Jdag, wmu, dJmax, dJ, *J20dag,
      musq, threemu1, threemu2, *J, *J20, *lambda, sign, lambda_gas, lambda_prv,
      lambda_nxt, fac, dl, frac;
  ActiveSet *as;
  CMO_PROF_FUNC_START();
  CMO_PROF_REGION_START("ALLOC & BOOLS");

  /* --- Retrieve active set as of transitions at wavelength nspect - */

  as = &spectrum.as[nspect];
  alloc_as(nspect, eval_operator, threadId);

  /* --- Check whether current active set includes a bound-bound
         and/or polarized transition and/or angle-dependent PRD
         transition, and/or polarization through background scattering.
         Otherwise, only angle-independent opacity and source functions
         are needed --                                 -------------- */

  /* --- Check for bound-bound transition in active set -- ---------- */

  boundbound = containsBoundBound(as);

  /* --- Check for line with angle-dependent PRD in set -- ---------- */

  PRD_angle_dep =
      (containsPRDline(as) && input.PRD_angle_dep != PRD_ANGLE_INDEP);

  /* --- Allocate temporary space --                   -------------- */

  if (eval_operator)
    Psi = (double *)malloc(Nspace * sizeof(double));
  else
    Psi = NULL;

  I = (double *)malloc(Nspace * sizeof(double));
  S = (double *)malloc(Nspace * sizeof(double));
  chi = (double *)malloc(Nspace * sizeof(double));

  /* --- Store current mean intensity, initialize new one to zero - - */

  Jdag = (double *)malloc(Nspace * sizeof(double));
  J = spectrum.J[nspect];
  for (k = 0; k < Nspace; k++)
    Jdag[k] = J[k];

  if (spectrum.updateJ)
    for (k = 0; k < Nspace; k++)
      J[k] = 0.0;
  CMO_PROF_REGION_END("ALLOC & BOOLS");
  /* --- Case of angle-dependent opacity and source function -- ----- */

  for (mu = 0; mu < Nrays; mu++) {
    wmu = 0.5 * geometry.wmu[mu];

    for (to_obs = 0; to_obs <= 1; to_obs++) {
      initialize = (mu == 0 && to_obs == 0);

      CMO_PROF_REGION_START("BG");
      if (initialize || (atmos.backgrflags[nspect] & HAS_LINE)) {
        if (input.backgr_in_mem) {
          // loadBackground(nspect, mu, to_obs);
          cmo_load_background(nspect, mu, to_obs);
        } else {
          readBackground(nspect, mu, to_obs);
        }
      }
      CMO_PROF_REGION_END("BG");

      if (initialize || boundbound)
        Opacity(nspect, mu, to_obs, initialize, threadId);

      if (eval_operator)
        addtoCoupling(nspect, threadId);

      CMO_PROF_REGION_START("SourceFn Prep");
      for (k = 0; k < Nspace; k++) {
        chi[k] = as->chi[k] + as->chi_c[k];
        S[k] = as->eta[k] + as->eta_c[k] + as->sca_c[k] * Jdag[k];
      }

      for (k = 0; k < Nspace; k++)
        S[k] /= chi[k];
      CMO_PROF_REGION_END("SourceFn Prep");
      // if (input.S_interpolation == S_LINEAR) {
      //   // printf("Sinterpol=LINEAR\n");
        Piecewise_1D(nspect, mu, to_obs, chi, S, I, Psi);
      // } else if (input.S_interpolation == CUBIC_HERMITE) {
      //   // printf("Sinterpol=CUBIC_HERMITE\n");
      //   Piecewise_Hermite_1D(nspect, mu, to_obs, chi, S, I, Psi);
      // } else {
      //   PieceBezier_1D(nspect, mu, to_obs, chi, S, I, Psi);
      // }

      CMO_PROF_REGION_START("EVAL OP");
      if (eval_operator) {
        for (k = 0; k < Nspace; k++)
          Psi[k] /= chi[k];
        addtoGamma(nspect, wmu, I, Psi, threadId);
      }
      CMO_PROF_REGION_END("EVAL OP");

      if (spectrum.updateJ) {
        CMO_PROF_REGION_START("UPDATE J");

        /* --- Accumulate mean intensity and rates -- ----------- */

        for (k = 0; k < Nspace; k++)
          J[k] += wmu * I[k];
        addtoRates(nspect, mu, to_obs, wmu, I, redistribute, threadId);

        /* --- Accumulate gas-frame mean intensity ------------- */

        if (atmos.NPRDactive > 0 && input.PRD_angle_dep == PRD_ANGLE_APPROX &&
            atmos.Nrays > 1) {
          CMO_PROF_REGION_START("Reg: PRDH");
          for (k = 0; k < Nspace; k++) {

            lamuk = nspect * (atmos.Nrays * 2 * Nspace) +
                    mu * (2 * Nspace) + to_obs * (Nspace) + k;

            idx0 = (lamuk == 0) ? 0 : spectrum.nc[lamuk - 1];

            for (idx = idx0; idx < spectrum.nc[lamuk]; idx++)
              spectrum.Jgas[spectrum.iprdh[idx]][k] +=
                  wmu * spectrum.cprdh[idx] * I[k];
          }
          CMO_PROF_REGION_END("Reg: PRDH");
        }   // Jgas accumulation endif

        if (containsPRDline(as) && input.PRD_angle_dep == PRD_ANGLE_DEP)
          writeImu(nspect, mu, to_obs, I);
        CMO_PROF_REGION_END("UPDATE J");
      }
      // NOTE(cmo): Need to save everything from here
      // Right, so here's the deal. Knowing that these matrices exist (and also
      // the ones for JDepth, and PsiStarDepth), if we change the running order
      // so that they are linear in k, we can simply use them to directly
      // compute the formal solution
      CMO_PROF_REGION_START("DEPTH DATA");
      for (int k = 0; k < Nspace; ++k)
      {
        int muk = k * atmos.Nrays * 2 + mu * 2 + to_obs;
        spectrum.IDepth[nspect][muk] = I[k];
        spectrum.SDepth[nspect][muk] = S[k];
        spectrum.chiDepth[nspect][muk] = chi[k];
      }
      CMO_PROF_REGION_END("DEPTH DATA");
    }

    /* --- Save emergent intensity --              -------------- */

    spectrum.I[nspect][mu] = I[0];
  }

  /* --- Write new J for current position in the spectrum -- -------- */

  CMO_PROF_REGION_START("DJ");
  dJmax = 0.0;
  if (spectrum.updateJ) {
    for (k = 0; k < Nspace; k++) {
      dJ = fabs(1.0 - Jdag[k] / J[k]);
      dJmax = MAX(dJmax, dJ);
    }
  }
  CMO_PROF_REGION_END("DJ");
  /* --- Clean up --                                 ---------------- */

  CMO_PROF_REGION_START("Reg: CleanUp");
  free_as(nspect, eval_operator, threadId);
  if (eval_operator)
    free(Psi);

  free(chi);
  free(I);
  free(S);

  free(Jdag);

  CMO_PROF_REGION_END("Reg: CleanUp");
  CMO_PROF_FUNC_END();
  return dJmax;
}
/* ------- begin -------------------------- Formal.c ---------------- */

double Formal(int nspect, bool_t eval_operator, bool_t redistribute, int threadId) {
  const char routineName[] = "Formal";
  register int k, mu, n;

  bool_t initialize, boundbound, polarized_as, polarized_c, PRD_angle_dep,
      to_obs, solveStokes, angle_dep;
  enum FeautrierOrder F_order;
  int Nrays = atmos.Nrays, lamuk;
  long Nspace = atmos.Nspace;
  long int idx, idx0;
  double *I, *chi, *S, **Ipol, **Spol, *Psi, *Jdag, wmu, dJmax, dJ, *J20dag,
      musq, threemu1, threemu2, *J, *J20, *lambda, sign, lambda_gas, lambda_prv,
      lambda_nxt, fac, dl, frac;
  ActiveSet *as;
  CMO_PROF_FUNC_START();

  /* --- Retrieve active set as of transitions at wavelength nspect - */
  CMO_PROF_REGION_START("ALLOC & BOOLS");
  as = &spectrum.as[nspect];
  alloc_as(nspect, eval_operator, threadId);

  /* --- Check whether current active set includes a bound-bound
         and/or polarized transition and/or angle-dependent PRD
         transition, and/or polarization through background scattering.
         Otherwise, only angle-independent opacity and source functions
         are needed --                                 -------------- */

  /* --- Check for bound-bound transition in active set -- ---------- */

  boundbound = containsBoundBound(as);

  /* --- Check for line with angle-dependent PRD in set -- ---------- */

  PRD_angle_dep =
      (containsPRDline(as) && input.PRD_angle_dep != PRD_ANGLE_INDEP);

  /* --- Check for polarized bound-bound transition in active set - - */

  polarized_as = containsPolarized(as);

  /* --- Check for polarized bound-bound transition in background - - */

  // polarized_c = atmos.backgrflags[nspect].ispolarized;
  polarized_c = atmos.backgrflags[nspect] & IS_POLARIZED;

  /* --- Determine if we solve for I, or for I, Q, U, V -- ---------- */

  solveStokes = (input.StokesMode == FULL_STOKES &&
                 (polarized_as || polarized_c || input.backgr_pol));

  /* --- Determine if we have to do angle-dependent opacity and
         emissivity --                                 -------------- */

  angle_dep =
      (polarized_as || polarized_c || PRD_angle_dep != PRD_ANGLE_INDEP ||
       (input.backgr_pol && input.StokesMode == FULL_STOKES) ||
      //  (atmos.moving && (boundbound || atmos.backgrflags[nspect].hasline)));
       (atmos.moving && (boundbound || (atmos.backgrflags[nspect] & HAS_LINE))));

  /* --- Allocate temporary space --                   -------------- */

  if (eval_operator)
    Psi = (double *)malloc(Nspace * sizeof(double));
  else
    Psi = NULL;

  if (solveStokes) {
    Ipol = matrix_double(4, Nspace);
    I = Ipol[0];
    Spol = matrix_double(4, Nspace);
    S = Spol[0];
  } else {
    I = (double *)malloc(Nspace * sizeof(double));
    S = (double *)malloc(Nspace * sizeof(double));
  }
  chi = (double *)malloc(Nspace * sizeof(double));

  /* --- Store current mean intensity, initialize new one to zero - - */

  Jdag = (double *)malloc(Nspace * sizeof(double));
  if (input.limit_memory) {
    J = (double *)malloc(Nspace * sizeof(double));
    readJlambda(nspect, Jdag);
  } else {
    J = spectrum.J[nspect];
    for (k = 0; k < Nspace; k++)
      Jdag[k] = J[k];
  }

  if (spectrum.updateJ)
    for (k = 0; k < Nspace; k++)
      J[k] = 0.0;
  CMO_PROF_REGION_END("ALLOC & BOOLS");

  /* --- Store current anisotropy, initialize new one to zero ---- -- */

  if (input.backgr_pol) {
    J20dag = (double *)malloc(Nspace * sizeof(double));
    if (input.limit_memory) {
      J20 = (double *)malloc(Nspace * sizeof(double));
      readJ20lambda(nspect, J20dag);
    } else {
      J20 = spectrum.J20[nspect];
      for (k = 0; k < Nspace; k++)
        J20dag[k] = J20[k];
    }
    if (spectrum.updateJ)
      for (k = 0; k < Nspace; k++)
        J20[k] = 0.0;
  }
  /* --- Case of angle-dependent opacity and source function -- ----- */

  if (angle_dep) {
    CMO_PROF_REGION_START("ANGLE DEP");
    for (mu = 0; mu < Nrays; mu++) {
      wmu = 0.5 * geometry.wmu[mu];
      if (input.backgr_pol) {
        musq = SQ(geometry.muz[mu]);
        threemu1 = TWOSQRTTWO * (3.0 * musq - 1.0);
        threemu2 = (3.0 * TWOSQRTTWO) * (musq - 1.0);
      }
      for (to_obs = 0; to_obs <= 1; to_obs++) {
        initialize = (mu == 0 && to_obs == 0);

        CMO_PROF_REGION_START("BG");
        if (initialize || (atmos.backgrflags[nspect] & HAS_LINE)) {
          if (input.backgr_in_mem) {
            // loadBackground(nspect, mu, to_obs);
            cmo_load_background(nspect, mu, to_obs);
          } else {
            readBackground(nspect, mu, to_obs);
          }
        }
        CMO_PROF_REGION_END("BG");

        if (initialize || boundbound)
          Opacity(nspect, mu, to_obs, initialize, threadId);

        CMO_PROF_REGION_START("SourceFn Prep");
        if (eval_operator)
          addtoCoupling(nspect, threadId);
        for (k = 0; k < Nspace; k++) {
          chi[k] = as->chi[k] + as->chi_c[k];
          S[k] = as->eta[k] + as->eta_c[k] + as->sca_c[k] * Jdag[k];
        }

        CMO_PROF_REGION_END("SourceFn Prep");

        if (solveStokes) {
          for (k = Nspace; k < 4 * Nspace; k++)
            Spol[0][k] = 0.0;

          /* --- Add emissivity due to active set for Q, U, V -- ---- */

          if (polarized_as) {
            for (k = Nspace; k < 4 * Nspace; k++)
              Spol[0][k] += as->eta[k];
          }
          /* --- Add emissivity due to background lines -- ---------- */

          if (polarized_c) {
            for (k = Nspace; k < 4 * Nspace; k++)
              Spol[0][k] += as->eta_c[k];
          }
          /* --- Add emissivity due to background scattering -- ----- */

          if (input.backgr_pol && input.StokesMode == FULL_STOKES) {
            for (k = 0; k < Nspace; k++) {
              Spol[0][k] += threemu1 * as->sca_c[k] * J20dag[k];
              Spol[1][k] += threemu2 * as->sca_c[k] * J20dag[k];
            }
          }
          for (n = 0; n < 4; n++) {
            for (k = 0; k < Nspace; k++)
              Spol[n][k] /= chi[k];
          }
          PiecewiseStokes(nspect, mu, to_obs, chi, Spol, Ipol, Psi);

        } else {
          CMO_PROF_REGION_START("DIV");
          for (k = 0; k < Nspace; k++)
            S[k] /= chi[k];
          CMO_PROF_REGION_END("DIV");
          if (input.S_interpolation == S_LINEAR) {
            // printf("Sinterpol=LINEAR\n");
            Piecewise_1D(nspect, mu, to_obs, chi, S, I, Psi);
          } else if (input.S_interpolation == CUBIC_HERMITE) {
            // printf("Sinterpol=CUBIC_HERMITE\n");
            Piecewise_Hermite_1D(nspect, mu, to_obs, chi, S, I, Psi);
          } else {
            PieceBezier_1D(nspect, mu, to_obs, chi, S, I, Psi);
          }
        }

        if (eval_operator) {
          CMO_PROF_REGION_START("EVAL OP");
          for (k = 0; k < Nspace; k++)
            Psi[k] /= chi[k];
          addtoGamma(nspect, wmu, I, Psi, threadId);
          CMO_PROF_REGION_END("EVAL OP");
        }

        if (spectrum.updateJ) {
          CMO_PROF_REGION_START("UPDATE J");

          /* --- Accumulate mean intensity and rates -- ----------- */

          for (k = 0; k < Nspace; k++)
            J[k] += wmu * I[k];
          addtoRates(nspect, mu, to_obs, wmu, I, redistribute, threadId);

          /* --- Accumulate anisotropy --            -------------- */

          if (input.backgr_pol) {
            for (k = 0; k < Nspace; k++)
              J20[k] += (threemu1 * Ipol[0][k] + threemu2 * Ipol[1][k]) * wmu;
          }

          /* --- Accumulate gas-frame mean intensity ------------- */

          if (atmos.NPRDactive > 0 && input.PRD_angle_dep == PRD_ANGLE_APPROX &&
              atmos.Nrays > 1) {
            CMO_PROF_REGION_START("Reg: PRDH");
            if (input.prdh_limit_mem) {

              sign = (to_obs) ? 1.0 : -1.0;

              for (k = 0; k < Nspace; k++) {

                // Observer's frame wavelenght grid
                lambda = spectrum.lambda;

                // previous, current and next wavelength shifted to gas rest
                // frame
                fac = (1. + spectrum.v_los[mu][k] * sign / CLIGHT);
                lambda_prv = lambda[MAX(nspect - 1, 0)] * fac;
                lambda_gas = lambda[nspect] * fac;
                lambda_nxt = lambda[MIN(nspect + 1, spectrum.Nspect - 1)] * fac;

                // do lambda_prv and lambda_gas bracket lambda points?
                if (lambda_prv != lambda_gas) {

                  dl = lambda_gas - lambda_prv;
                  for (idx = 0; idx < spectrum.Nspect; idx++) {
                    if (lambda[idx] > lambda_prv && lambda[idx] <= lambda_gas) {
                      frac = (lambda[idx] - lambda_prv) / dl;
                      spectrum.Jgas[idx][k] += frac * wmu * I[k];
                    }
                  }

                } else {

                  // edge case, use constant extrapolation for
                  // lambda[idx]<lambda gas
                  for (idx = 0; idx < spectrum.Nspect; idx++) {
                    if (lambda[idx] < lambda_gas)
                      spectrum.Jgas[idx][k] += wmu * I[k];
                  }
                }

                // do lambda_gas and lambda_nxt bracket lambda points?
                if (lambda_gas != lambda_nxt) {

                  dl = lambda_nxt - lambda_gas;
                  for (idx = 0; idx < spectrum.Nspect; idx++) {
                    if (lambda[idx] > lambda_gas && lambda[idx] < lambda_nxt) {
                      frac = (lambda[idx] - lambda_gas) / dl;
                      spectrum.Jgas[idx][k] += (1.0 - frac) * wmu * I[k];
                    }
                  }

                } else {
                  // edge case, use constant extrapolation for
                  // lambda[idx]>lambda gas
                  for (idx = 0; idx < spectrum.Nspect; idx++) {
                    if (lambda[idx] > lambda_gas)
                      spectrum.Jgas[idx][k] += wmu * I[k];
                  }
                }

              } // spatial location

            } else {

              for (k = 0; k < Nspace; k++) {

                lamuk = nspect * (atmos.Nrays * 2 * Nspace) +
                        mu * (2 * Nspace) + to_obs * (Nspace) + k;

                idx0 = (lamuk == 0) ? 0 : spectrum.nc[lamuk - 1];

                for (idx = idx0; idx < spectrum.nc[lamuk]; idx++)
                  spectrum.Jgas[spectrum.iprdh[idx]][k] +=
                      wmu * spectrum.cprdh[idx] * I[k];
              }

            } // prdh_limit_mem switch
            CMO_PROF_REGION_END("Reg: PRDH");
          }   // Jgas accumulation endif

          if (containsPRDline(as) && input.PRD_angle_dep == PRD_ANGLE_DEP)
            writeImu(nspect, mu, to_obs, I);
          CMO_PROF_REGION_END("UPDATE J");
        }
        // NOTE(cmo): Need to save everything from here
        CMO_PROF_REGION_START("DEPTH DATA");
        for (int k = 0; k < Nspace; ++k)
        {
          int muk = k * atmos.Nrays * 2 + mu * 2 + to_obs;
          spectrum.IDepth[nspect][muk] = I[k];
          spectrum.SDepth[nspect][muk] = S[k];
          spectrum.chiDepth[nspect][muk] = chi[k];
        }
        CMO_PROF_REGION_END("DEPTH DATA");
      }

      /* --- Save emergent intensity --              -------------- */

      spectrum.I[nspect][mu] = I[0];
      if (solveStokes) {
        spectrum.Stokes_Q[nspect][mu] = Ipol[1][0];
        spectrum.Stokes_U[nspect][mu] = Ipol[2][0];
        spectrum.Stokes_V[nspect][mu] = Ipol[3][0];
      }
    }

    CMO_PROF_REGION_END("ANGLE DEP");
  } else {

    CMO_PROF_REGION_START("ANGLE INDEP");
    /* --- The angle-independent case --               -------------- */
    CMO_PROF_REGION_START("BG");
    if (input.backgr_in_mem) {
      // loadBackground(nspect, 0, 0);
      cmo_load_background(nspect, 0, 0);
    } else {
      readBackground(nspect, 0, 0);
    }
    CMO_PROF_REGION_END("BG");

    Opacity(nspect, 0, 0, initialize = TRUE, threadId);
    if (eval_operator)
      addtoCoupling(nspect, threadId);

    CMO_PROF_REGION_START("SourceFn Prep");
    for (k = 0; k < Nspace; k++) {
      chi[k] = as->chi[k] + as->chi_c[k];
      S[k] = (as->eta[k] + as->eta_c[k] + as->sca_c[k] * Jdag[k]) / chi[k];
    }
    CMO_PROF_REGION_END("SourceFn Prep");

    for (mu = 0; mu < Nrays; mu++) {
      spectrum.I[nspect][mu] =
          Feautrier(nspect, mu, chi, S, F_order = STANDARD, I, Psi);
      if (eval_operator) {
        CMO_PROF_REGION_START("DIV");
        for (k = 0; k < Nspace; k++)
          Psi[k] /= chi[k];
        CMO_PROF_REGION_END("DIV");
        addtoGamma(nspect, geometry.wmu[mu], I, Psi, threadId);
      }

      if (spectrum.updateJ) {
        CMO_PROF_REGION_START("UPDATE J INDEP");
        for (k = 0; k < Nspace; k++)
          J[k] += I[k] * geometry.wmu[mu];
        addtoRates(nspect, mu, 0, geometry.wmu[mu], I, redistribute, threadId);

        /* --- Accumulate gas-frame mean intensity, which is the same
           as J in the angle-independent case ------------- */
        if (atmos.NPRDactive > 0 && input.PRD_angle_dep == PRD_ANGLE_APPROX) {
          for (k = 0; k < Nspace; k++)
            spectrum.Jgas[nspect][k] += I[k] * geometry.wmu[mu];
        }
        CMO_PROF_REGION_END("UPDATE J INDEP");
      }
      // NOTE(cmo): The I computed in the Feautrier solution is actually 
      // P = 0.5 * (I+ + I-), so we can assign it to both to_obs directions.
      CMO_PROF_REGION_START("DEPTH DATA");
      for (int to_obs = 0; to_obs <= 1; ++to_obs)
      {
        for (int k = 0; k < Nspace; ++k)
        {
          int muk = k * atmos.Nrays * 2 + mu * 2 + to_obs;
          spectrum.IDepth[nspect][muk] = I[k];
          spectrum.SDepth[nspect][muk] = S[k];
          spectrum.chiDepth[nspect][muk] = chi[k];
        }
      }
      CMO_PROF_REGION_END("DEPTH DATA");
    }

    CMO_PROF_REGION_END("ANGLE INDEP");
  }

  /* --- Write new J for current position in the spectrum -- -------- */

  CMO_PROF_REGION_START("DJ");
  dJmax = 0.0;
  if (spectrum.updateJ) {
    for (k = 0; k < Nspace; k++) {
      dJ = fabs(1.0 - Jdag[k] / J[k]);
      dJmax = MAX(dJmax, dJ);
    }
    if (input.limit_memory) {
      writeJlambda(nspect, J);
      if (input.backgr_pol)
        writeJ20lambda(nspect, J20);
    }
  }
  CMO_PROF_REGION_END("DJ");
  /* --- Clean up --                                 ---------------- */

  CMO_PROF_REGION_START("Reg: CleanUp");
  free_as(nspect, eval_operator, threadId);
  if (eval_operator)
    free(Psi);

  free(chi);
  if (solveStokes) {
    freeMatrix((void **)Ipol);
    freeMatrix((void **)Spol);
  } else {
    free(I);
    free(S);
  }

  free(Jdag);
  if (input.limit_memory)
    free(J);
  if (input.backgr_pol) {
    free(J20dag);
    if (input.limit_memory)
      free(J20);
  }

  CMO_PROF_REGION_END("Reg: CleanUp");
  CMO_PROF_FUNC_END();
  return dJmax;
}
/* ------- end ---------------------------- Formal.c ---------------- */
