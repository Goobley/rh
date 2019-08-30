/* ------- file: -------------------------- background.c ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Jul 24 12:52:46 2013 --

       --------------------------                      ----------RH-- */

/* Driving subroutine for background opacity sources.

 * Included at the moment:

  ++ Thomson scattering by free electrons

  ++ Hydrogen:
    -- Bound-free absorption and emission by neutral Hydrogen
    -- Free-free absorption and emission by neutral Hydrogen
    -- Rayleigh scattering by neutral Hydrogen and Helium
    -- Rayleigh scattering by molecular Hydrogen (H2)
    -- Bound-free absorption and emission by Hminus (H^-)
    -- Free-free absorption and emission by Hminus (H + e)
    -- Free-free absorption and emission by H2minus (H2 + e)
    -- Free-free absorption and emission by H2plus (H + H^+)

  ++ Metals:
    -- Bound-free absorption and emission from metals specified in
       file background.input.
    -- Bound-bound absorption and emission from metals specified in
       file background.input.
    -- LTE Bound-bound absorption and emission all elements from
       a Kurucz line list.

  ++ Molecules:
    -- Chemical equilibrium is calculated for molecules specified
       in file background.input and populations of constituent atoms
       are reduced accordingly.
    -- molecular opacities (LTE) may be taken into account by specifying
       data files with transition lists in the molecular input files.
    -- Bound-free absorption and emission by OH and CH molecules.

 * Atomic models are specified in atoms.input, molecules in
   molecules.input

 * Entries for the atoms.input and molecules.input files should have
   the form, respectively:

    -------------------------------------------------------------------
   |                                                                  |
   |   Nmetal                                                         |
   |                                                                  |
   |   model file ACTIVE/PASSIVE  INITIAL_SOLUTION   population file  |
   |                          .                                       |
   |                          .                                       |
   |                                                                  |
   |   Nmolecule                                                      |
   |                                                                  |
   |   model file ACTIVE/PASSIVE  INITIAL_SOLUTION                    |
   |                 .                                                |
   |                 .                                                |
    -------------------------------------------------------------------


   Nmetal and Nmolecule are the number of metal and molecule entries.
   metalID is the two-character atomID, the next entry is either
   set to LTE or NLTE, model_file is the input file containing atomic
   data for this metal (generic atomic input data format), and
   population_file is the input file containing the NLTE population
   numbers from a previous calculation. This last entry is only read when
   the second entry is set to NLTE.

   -- Units:
      Wavelengths are given in nm, densities in m^-3, opacities in m^2,
      and emissivities in J s^-1 Hz^-1 sr^-1.

 Note: The model atom file for hydrogen is specified in keyword.input.
       If H_LTE = TRUE is specified there LTE hydrogen populations are
       used. See: distribute_nH in the file hydrogen.c

 Note: Scattering opacity is added to total opacity after all
       contributions have been computed.

 Note: The quantities chi_ai and eta_ai store the angle-inpendent
       opacities and emissivities in case atmos.moving == TRUE.
       In static atmospheres these quantities are just mapped to
       atmos.chi_c and atmos.eta_c to save memory space.

 Note: If write_analyze_output == FALSE the auxiliary output files for
       the Background Record Structure (BRS), metals, and molecules
       are NOT written. This option is used when Background is called
       from programs like solveray (formal solution along one specific
       ray) in cases with moving atmospheres (angle-dependent opacity).

 Note: If equilibria_only is set to TRUE only the electron density,
       LTE populations and collisions, and chemical equilibria are
       evaluated.

 Note: Record numbers stored in atmos.backgrrecno refer to records
       of the size atmos.Nspace. If a polarized line is present 9
       (4 + 4 + 1, no magneto-optical effects), or 12 (7 + 4 +1, with
       magneto-optical effects) records are used, otherwise 3 (1 + 1 + 1)
       records.
       If the atmosphere is moving, or if a polarized line is present
       data is stored for each angle and wavelength, otherwise data
       is stored once for each wavelength only.
       --                                              -------------- */

#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "constant.h"
#include "background.h"
#include "error.h"
#include "statistics.h"
#include "inputs.h"
#include "xdr.h"
#include "CmoProfile.h"

#define COMMENT_CHAR "#"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];

void cmo_init_background()
{
  atmos.Nrlk = 0;
  readKuruczLines(input.KuruczData);
  if (atmos.Nrlk > 0) {
    qsort(atmos.rlk_lines, atmos.Nrlk, sizeof(RLK_Line), rlk_ascend);
  }
  atmos.backgrflags = (flags *)malloc(spectrum.Nspect * sizeof(flags));
  if (atmos.moving || atmos.Stokes) {
    atmos.backgrrecno =
        (long *)calloc(2 * spectrum.Nspect * atmos.Nrays, sizeof(long));
  } else
    atmos.backgrrecno = (long *)calloc(spectrum.Nspect, sizeof(long));
}

typedef struct BackgroundState
{
  double* chi_ai;
  double* eta_ai;
  double* sca_ai;

  double* chi_c;
  double* eta_c;
  double* sca_c;
  double* chip;
  double* chip_c;

  double* chi;
  double* eta;
  double* scatt;

  double* Bnu;

  Atom* He;
  int NrecStokes;
  bool_t do_fudge;
  int Nfudge;
  double* lambda_fudge;
  double** fudge;
  double* thomson;

  PassiveBbState* bbState;

} BackgroundState;

static void inner_background(BackgroundState* bg, int nspect, int threadId)
{
  // input.backgr_in_mem = TRUE;

  double wavelength = spectrum.lambda[nspect];
  // double Bnu[atmos.Nspace];

  /* --- The Planck function at this wavelength --   -------------- */

  Planck(atmos.Nspace, atmos.T, wavelength, bg->Bnu);

  /* --- Initialize the flags for this wavelength -- -------------- */

  // atmos.backgrflags[nspect].hasline = FALSE;
  // atmos.backgrflags[nspect].ispolarized = FALSE;
  atmos.backgrflags[nspect] = 0;

  /* --- Initialize angle-independent quantities --  -------------- */

  for (int k = 0; k < atmos.Nspace; k++) {
    bg->chi_ai[k] = 0.0;
    bg->eta_ai[k] = 0.0;
    bg->sca_ai[k] = bg->thomson[k];
  }
  /* --- Negative hydrogen ion, bound-free and free-free -- ------- */

  if (cmo_Hminus_bf(&input.splineState[threadId], wavelength, bg->chi, bg->eta)) {
    for (int k = 0; k < atmos.Nspace; k++) {
      bg->chi_ai[k] += bg->chi[k];
      bg->eta_ai[k] += bg->eta[k];
    }
  }
  if (cmo_Hminus_ff(wavelength, bg->chi)) {
    for (int k = 0; k < atmos.Nspace; k++) {
      bg->chi_ai[k] += bg->chi[k];
      bg->eta_ai[k] += bg->chi[k] * bg->Bnu[k];
    }
  }
  /* --- Opacity fudge factors, applied to Hminus opacity -- ------ */

  if (bg->do_fudge) {
    double Hmin_fudge;
    Linear(bg->Nfudge, bg->lambda_fudge, bg->fudge[0], 
           1, &wavelength, &Hmin_fudge, FALSE);
    for (int k = 0; k < atmos.Nspace; k++) {
      bg->chi_ai[k] *= Hmin_fudge;
      bg->eta_ai[k] *= Hmin_fudge;
    }
  }
  /* --- Opacities from bound-free transitions in OH and CH -- ---- */

  // if (FALSE)
  // {
  if (OH_bf_opac(wavelength, bg->chi, bg->eta)) {
    for (int k = 0; k < atmos.Nspace; k++) {
      bg->chi_ai[k] += bg->chi[k];
      bg->eta_ai[k] += bg->eta[k];
    }
  }
  if (CH_bf_opac(wavelength, bg->chi, bg->eta)) {
    for (int k = 0; k < atmos.Nspace; k++) {
      bg->chi_ai[k] += bg->chi[k];
      bg->eta_ai[k] += bg->eta[k];
    }
  }
  /* --- Neutral Hydrogen Bound-Free and Free-Free --  ------------ */

  if (Hydrogen_bf(wavelength, bg->chi, bg->eta)) {
    for (int k = 0; k < atmos.Nspace; k++) {
      bg->chi_ai[k] += bg->chi[k];
      bg->eta_ai[k] += bg->eta[k];
    }
  }
  Hydrogen_ff(wavelength, bg->chi);
  for (int k = 0; k < atmos.Nspace; k++) {
    bg->chi_ai[k] += bg->chi[k];
    bg->eta_ai[k] += bg->chi[k] * bg->Bnu[k];
  }
  /* --- Rayleigh scattering by neutral hydrogen --  -------------- */

  if (Rayleigh(wavelength, atmos.H, bg->scatt)) {
    for (int k = 0; k < atmos.Nspace; k++) {
      bg->sca_ai[k] += bg->scatt[k];
    }
  }
  /* --- Rayleigh scattering by neutral helium --    -------------- */

  if (bg->He && Rayleigh(wavelength, bg->He, bg->scatt)) {
    for (int k = 0; k < atmos.Nspace; k++) {
      bg->sca_ai[k] += bg->scatt[k];
    }
  }
  /* --- Absorption by H + H^+ (referred to as H2plus free-free) -- */

  if (cmo_H2plus_ff(wavelength, bg->chi)) {
    for (int k = 0; k < atmos.Nspace; k++) {
      bg->chi_ai[k] += bg->chi[k];
      bg->eta_ai[k] += bg->chi[k] * bg->Bnu[k];
    }
  }
  /* --- Rayleigh scattering and free-free absorption by
          molecular hydrogen --                       -------------- */

  // }
  if (Rayleigh_H2(wavelength, bg->scatt)) {
    for (int k = 0; k < atmos.Nspace; k++) {
      bg->sca_ai[k] += bg->scatt[k];
    }
  }
  // if (FALSE)
  // {
  if (cmo_H2minus_ff(wavelength, bg->chi)) {
    for (int k = 0; k < atmos.Nspace; k++) {
      bg->chi_ai[k] += bg->chi[k];
      bg->eta_ai[k] += bg->chi[k] * bg->Bnu[k];
    }
  }
  // }
  /* --- Bound-Free opacities due to ``metals'' --   -------------- */

  double metal_fudge = 1.0;
  if (bg->do_fudge) {
    Linear(bg->Nfudge, bg->lambda_fudge, bg->fudge[2], 
           1, &wavelength, &metal_fudge, FALSE);
  }
  /* --- Note: Hydrogen bound-free opacities are calculated in
          routine Hydrogen_bf --                      -------------- */

  // Metal_bf(wavelength, atmos.Natom - 1, atmos.atoms + 1, bg->chi, bg->eta);
  cmo_Metal_bf(&input.splineState[threadId], wavelength, 
               atmos.Natom - 1, atmos.atoms + 1, bg->chi, bg->eta);
  for (int k = 0; k < atmos.Nspace; k++) {
    bg->chi_ai[k] += bg->chi[k] * metal_fudge;
    bg->eta_ai[k] += bg->eta[k] * metal_fudge;
  }
  /* --- Add the scattering opacity to the absorption part to store
          the total opacity --                        -------------- */

  double scatt_fudge = 1.0;
  if (bg->do_fudge) {
    Linear(bg->Nfudge, bg->lambda_fudge, bg->fudge[1], 
           1, &wavelength, &scatt_fudge, FALSE);
  }

  for (int k = 0; k < atmos.Nspace; k++) {
    bg->sca_ai[k] *= scatt_fudge;
    bg->chi_ai[k] += bg->sca_ai[k];
  }
  /* --- Now the contributions that may be angle-dependent due to the
          presence of atomic or molecular lines --    -------------- */

  if (atmos.moving || atmos.Stokes) {
    for (int mu = 0; mu < atmos.Nrays; mu++) {
      for (int to_obs = 0; to_obs <= 1; to_obs++) {
        int index = 2 * (nspect * atmos.Nrays + mu) + to_obs;

        /* --- First, copy the angle-independent parts -- --------- */

        for (int k = 0; k < atmos.Nspace; k++) {
          bg->chi_c[k] = bg->chi_ai[k];
          bg->eta_c[k] = bg->eta_ai[k];
          bg->sca_c[k] = bg->sca_ai[k];
        }

        /* --- Zero the polarized quantities, if necessary -- ----- */

        if (atmos.Stokes) {
          for (int k = atmos.Nspace; k < 4 * atmos.Nspace; k++) {
            bg->chi_c[k] = 0.0;
            bg->eta_c[k] = 0.0;
          }
          if (input.magneto_optical)
            for (int k = 0; k < 3 * atmos.Nspace; k++)
              bg->chip_c[k] = 0.0;
        }
        /* --- Add opacity from passive atomic lines (including
                hydrogen) --                          -------------- */

        if (input.allow_passive_bb) {
          flags backgrflags =
              passive_bb_stateless(bg->bbState, wavelength, nspect, mu, to_obs, bg->chi, bg->eta, bg->chip);
              // passive_bb(wavelength, nspect, mu, to_obs, bg->chi, bg->eta, bg->chip);

          int nStokes = 1;
          // if (backgrflags.hasline) {
          if (backgrflags & HAS_LINE) {
            // atmos.backgrflags[nspect].hasline = TRUE;
            atmos.backgrflags[nspect] |= HAS_LINE;
            if (backgrflags & IS_POLARIZED) {
              nStokes = 4;
              atmos.backgrflags[nspect] |= IS_POLARIZED;
              if (input.magneto_optical) {
                for (int k = 0; k < 3 * atmos.Nspace; k++)
                  bg->chip_c[k] += bg->chip[k];
              }
            } 

            for (int k = 0; k < nStokes * atmos.Nspace; k++) {
              bg->chi_c[k] += bg->chi[k];
              bg->eta_c[k] += bg->eta[k];
            }
          }
        }
        /* --- Add opacity from Kurucz line list --  -------------- */

        if (atmos.Nrlk > 0) {
          flags backgrflags = rlk_opacity(wavelength, nspect, mu, to_obs, 
                                          bg->chi, bg->eta, bg->scatt, bg->chip);
          int nStokes = 1;
          if (backgrflags & HAS_LINE) {
            // atmos.backgrflags[nspect].hasline = TRUE;
            atmos.backgrflags[nspect] |= HAS_LINE;
            // if (backgrflags.ispolarized) {
            if (backgrflags & IS_POLARIZED) {
              nStokes = 4;
              // atmos.backgrflags[nspect].ispolarized = TRUE;
              atmos.backgrflags[nspect] |= IS_POLARIZED;
              if (input.magneto_optical) {
                for (int k = 0; k < 3 * atmos.Nspace; k++)
                  bg->chip_c[k] += bg->chip[k];
              }
            }

            for (int k = 0; k < nStokes * atmos.Nspace; k++) {
              bg->chi_c[k] += bg->chi[k];
              bg->eta_c[k] += bg->eta[k];
            }
            if (input.rlkscatter) {
              for (int k = 0; k < atmos.Nspace; k++) {
                bg->sca_c[k] += bg->scatt[k];
                bg->chi_c[k] += bg->scatt[k];
              }
            }
          }
        }
        /* --- Add opacity from molecular lines --   -------------- */

        flags backgrflags =
            MolecularOpacity(wavelength, nspect, mu, to_obs, 
                             bg->chi, bg->eta, bg->chip);
        int nStokes = 1;
        if (backgrflags & HAS_LINE) {
          // atmos.backgrflags[nspect].hasline = TRUE;
          atmos.backgrflags[nspect] |= HAS_LINE;
          if (backgrflags & IS_POLARIZED) {
            nStokes = 4;
            // atmos.backgrflags[nspect].ispolarized = TRUE;
            atmos.backgrflags[nspect] |= IS_POLARIZED;
            if (input.magneto_optical) {
              for (int k = 0; k < 3 * atmos.Nspace; k++)
                bg->chip_c[k] += bg->chip[k];
            }
          }

          for (int k = 0; k < nStokes * atmos.Nspace; k++) {
            bg->chi_c[k] += bg->chi[k];
            bg->eta_c[k] += bg->eta[k];
          }
        }
        /* --- Store angle-dependent results only if at least one
                line was found at this wavelength --  -------------- */

        // This is very messy logic, but it works
        if ((mu == atmos.Nrays - 1 && to_obs) ||
            // (atmos.backgrflags[nspect].hasline &&
            ((atmos.backgrflags[nspect] & HAS_LINE) &&
              // (atmos.moving || atmos.backgrflags[nspect].ispolarized))) {
              (atmos.moving || (atmos.backgrflags[nspect] & IS_POLARIZED)))) {
          // If we've accumulated, or there was angle-dependent info

          // if ((mu == atmos.Nrays - 1 && to_obs) &&
          //     !((atmos.backgrflags[nspect] & HAS_LINE) &&
          //       (atmos.moving || (atmos.backgrflags[nspect] & IS_POLARIZED)))) {
          //   // If we've accumulated across all of the mus and there wasn't any angle-dependent info, store everything into mu=0, to_obs=0
          //   int lamu = nspect * (2 * atmos.Nrays);
          //   atmos.backgrrecno[lamu] += cmo_store_background(nspect, 0, 0,
          //                     bg->chi_c, bg->eta_c, bg->sca_c, bg->chip_c);
          // } else {
          //   // Otherwise put the angle-dependent info where it belongs
          //   int lamu = nspect * (2 * atmos.Nrays) + mu * 2 + to_obs;
          //   atmos.backgrrecno[lamu] += cmo_store_background(nspect, mu,
          //           to_obs, bg->chi_c, bg->eta_c, bg->sca_c, bg->chip_c);
          //     // mpi.backgrrecno += writeBackground(nspect, mu, to_obs, chi_c,
          //     //                                    eta_c, sca_c, chip_c);
          // }
          int lamu = nspect * (2 * atmos.Nrays) + mu * 2 + to_obs;
          atmos.backgrrecno[lamu] += cmo_store_background(nspect, mu,
                  to_obs, bg->chi_c, bg->eta_c, bg->sca_c, bg->chip_c);
        }
      }
    }
  } else {
    /* --- Angle-independent case. First, add opacity from passive
            atomic lines (including hydrogen) --      -------------- */
    if (input.allow_passive_bb) {
      flags backgrflags = passive_bb_stateless(bg->bbState, wavelength, nspect, 0, TRUE, 
                                     bg->chi, bg->eta, NULL);
      // if (backgrflags.hasline) {
      if (backgrflags & HAS_LINE) {
        atmos.backgrflags[nspect] |= HAS_LINE;
        for (int k = 0; k < atmos.Nspace; k++) {
          bg->chi_c[k] += bg->chi[k];
          bg->eta_c[k] += bg->eta[k];
        }
      }
    }

    /* --- Add opacity from Kurucz line list --      -------------- */
    if (atmos.Nrlk > 0) {
      flags backgrflags =
          rlk_opacity(wavelength, nspect, 0, TRUE, 
                      bg->chi, bg->eta, bg->scatt, NULL);
      // if (backgrflags.hasline) {
      if (backgrflags & HAS_LINE) {
        // atmos.backgrflags[nspect].hasline = TRUE;
        atmos.backgrflags[nspect] |= HAS_LINE;
        for (int k = 0; k < atmos.Nspace; k++) {
          bg->chi_c[k] += bg->chi[k];
          bg->eta_c[k] += bg->eta[k];
        }
        if (input.rlkscatter) {
          for (int k = 0; k < atmos.Nspace; k++) {
            bg->sca_c[k] += bg->scatt[k];
            bg->chi_c[k] += bg->scatt[k];
          }
        }
      }
    }
    /* --- Add opacity from molecular lines --       -------------- */

    flags backgrflags =
        MolecularOpacity(wavelength, nspect, 0, TRUE, 
                         bg->chi, bg->eta, NULL);
    // if (backgrflags.hasline) {
    if (backgrflags & HAS_LINE) {
      atmos.backgrflags[nspect] |= HAS_LINE;
      for (int k = 0; k < atmos.Nspace; k++) {
        bg->chi_c[k] += bg->chi[k];
        bg->eta_c[k] += bg->eta[k];
      }
    }
    /* --- Store results --                          -------------- */

    // atmos.backgrrecno[nspect] = backgrrecno;
    // backgrrecno += writeBackground(nspect, 0, 0, chi_c, eta_c, sca_c, NULL);
    atmos.backgrrecno[nspect] += cmo_store_background(nspect, 0, 0, 
                            bg->chi_c, bg->eta_c, bg->sca_c, NULL);
  }
}

static void inner_background_sched(void* userdata, struct scheduler* s, struct sched_task_partition range, sched_uint threadId)
{
  CMO_PROF_FUNC_START();
  BackgroundState* bg = &((BackgroundState*)userdata)[threadId];
  for (int nspect = range.start; nspect < range.end; ++nspect)
  {
    inner_background(bg, nspect, threadId);
  }
  CMO_PROF_FUNC_END();
}

void cmo_background(bool_t write_analyze_output, bool_t equilibria_only)
{
  CMO_PROF_FUNC_START();
  // input.backgr_in_mem = TRUE;

  static int ne_iter = 0;
  if (input.solve_ne == ONCE || input.solve_ne == ITERATION)
  {
    bool_t fromScratch = (input.solve_ne == ONCE ||
                   (input.solve_ne == ITERATION && ne_iter == 0))
                      ? TRUE
                      : FALSE;
    Solve_ne(atmos.ne, fromScratch);
    ++ne_iter;
  }
  SetLTEQuantities();

  if (input.NonICE)
  {
    readMolecules(MOLECULAR_CONCENTRATION_FILE);
  }
  else
  {
    ChemicalEquilibrium(N_MAX_CHEM_ITER, CHEM_ITER_LIMIT);
  }

  if (equilibria_only)
  {
    CMO_PROF_FUNC_END();
    return;
  }

  if (input.old_background) {
    if ((atmos.fd_background = open(input.background_File, O_RDONLY, 0)) ==
        -1) {
      sprintf(messageStr, "Unable to open input file %s",
              input.background_File);
      Error(ERROR_LEVEL_2, __func__, messageStr);
    }
    readBRS();
    CMO_PROF_FUNC_END();
    return;
  }

  bool_t do_fudge = FALSE;
  char inputLine[MAX_LINE_SIZE];
  int Nfudge;
  double* lambda_fudge;
  double** fudge;
  if (strcmp(input.fudgeData, "none")) {
    do_fudge = TRUE;
    bool_t exit_on_EOF;

    /* --- Read wavelength-dependent fudge factors to compensate for
           missing UV backround line haze --           -------------- */

    FILE* fp_fudge;
    if ((fp_fudge = fopen(input.fudgeData, "r")) == NULL) {
      sprintf(messageStr, "Unable to open input file %s", input.fudgeData);
      Error(ERROR_LEVEL_2, __func__, messageStr);
    }
    sprintf(messageStr, "\n-Fudging background opacities with file\n  %s\n\n",
            input.fudgeData);
    Error(MESSAGE, __func__, messageStr);

    getLine(fp_fudge, COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
    sscanf(inputLine, "%d", &Nfudge);
    lambda_fudge = (double *)malloc(Nfudge * sizeof(double));
    fudge = matrix_double(3, Nfudge);
    for (int n = 0; n < Nfudge; n++) {
      getLine(fp_fudge, COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
      sscanf(inputLine, "%lf %lf %lf %lf", &lambda_fudge[n], &fudge[0][n],
             &fudge[1][n], &fudge[2][n]);
    }
    for (int n = 0; n < 3 * Nfudge; n++)
      fudge[0][n] += 1.0;
    fclose(fp_fudge);
  }

  int NrecStokes = 1;
  if (atmos.Stokes)
    NrecStokes = 4;
  
  // if (!input.backgr_in_mem)
  // {
  //   assert(false && "CMO Background assumes backgr_in_mem");
  // }

  // TODO(cmo): Should really handle non-moving/non-Stokes case here
  int allocSize = spectrum.Nspect;
  if (atmos.moving || atmos.Stokes)
    allocSize *= 2 * atmos.Nrays;
  if (!atmos.chi_b)
  {
    atmos.chi_b = (double**)calloc(allocSize, sizeof(double*));
    atmos.eta_b = (double**)calloc(allocSize, sizeof(double*));
    atmos.sca_b = (double**)calloc(allocSize, sizeof(double*));
    if (input.magneto_optical)
    {
      atmos.chip_b = (double**)calloc(allocSize, sizeof(double*));
    }
  }

  if (atmos.chi_b[0])
  {
    if (atmos.moving || atmos.Stokes)
    {
      int prevIdx = 0;
      for (int i = 1; i < allocSize; ++i)
      {
        if (atmos.chi_b[i] == atmos.chi_b[prevIdx])
          continue;

        free(atmos.chi_b[prevIdx]);
        free(atmos.eta_b[prevIdx]);
        free(atmos.sca_b[prevIdx]);
        if (input.magneto_optical && atmos.chip_b[prevIdx])
          free(atmos.chip_b[prevIdx]);
        prevIdx = i;
      }
      free(atmos.chi_b[prevIdx]);
      free(atmos.eta_b[prevIdx]);
      free(atmos.sca_b[prevIdx]);
      if (input.magneto_optical && atmos.chip_b[prevIdx])
        free(atmos.chip_b[prevIdx]);
    }
    else
    {
      for (int i = 0; i < allocSize; ++i)
      {
        if (atmos.chi_b[i])
        {
          free(atmos.chi_b[i]);
          free(atmos.eta_b[i]);
          free(atmos.sca_b[i]);
          if (input.magneto_optical)
          {
            free(atmos.chip_b[i]);
          }
        }
      }
    }
    memset(atmos.chi_b, 0, allocSize * sizeof(double*));
    memset(atmos.eta_b, 0, allocSize * sizeof(double*));
    memset(atmos.sca_b, 0, allocSize * sizeof(double*));
    if (input.magneto_optical)
      memset(atmos.chip_b, 0, allocSize * sizeof(double*));
  }

  // Allocate these guys per thread
  /* --- Allocate temporary storage space. The quantities are used
         for the following purposes:

       - chi, eta, scatt: Get contributions to opacity, emissivity,
         and scattering opacity, respectively, from a specific process
         for a given wavelength and possibly angle.

       - chi_c, eta_c, sca_c: Collect the total opacity, emissivity
         and scattering opacity for a given wavelength and possibly
         angle.

       - chi_ai, eta_ai: Collect the angle-independent part of
         opacity and emissivity for each wavelength so that these
         need not be recalculated in an angle-dependent case.
         When the atmosphere is not moving and has no magnetic fields
         these just point to the total quantities chi_c and eta_c.

   Note: In case of magnetic fields in the atmosphere chi, eta and
         chip, and chi_c, eta_c and chip_c contain all four Stokes
         parameters, and should be allocated a 4 and 3 times larger
         storage space, respectively.
         --                                            -------------- */
  // double *chi_c, *eta_c, *sca_c, *chi, *eta, *scatt, *chip_c, *chip, 
  //        *chi_ai, *eta_ai, *sca_ai;
  // chi_c = (double *)malloc(NrecStokes * atmos.Nspace * sizeof(double));
  // eta_c = (double *)malloc(NrecStokes * atmos.Nspace * sizeof(double));
  // sca_c = (double *)malloc(atmos.Nspace * sizeof(double));

  // chi = (double *)malloc(NrecStokes * atmos.Nspace * sizeof(double));
  // eta = (double *)malloc(NrecStokes * atmos.Nspace * sizeof(double));
  // scatt = (double *)malloc(atmos.Nspace * sizeof(double));

  // if (atmos.Stokes && input.magneto_optical) {
  //   chip = (double *)malloc(3 * atmos.Nspace * sizeof(double));
  //   chip_c = (double *)malloc(3 * atmos.Nspace * sizeof(double));
  // } else {
  //   chip = NULL;
  //   chip_c = NULL;
  // }

  // if (atmos.moving || atmos.Stokes) {
  //   chi_ai = (double *)malloc(atmos.Nspace * sizeof(double));
  //   eta_ai = (double *)malloc(atmos.Nspace * sizeof(double));
  //   sca_ai = (double *)malloc(atmos.Nspace * sizeof(double));
  // } else {
  //   chi_ai = chi_c;
  //   eta_ai = eta_c;
  //   sca_ai = sca_c;
  // }
  /* --- Thomson scattering by free electrons is wavelength independent
         in non-relativistic limit so we compute it only once -- ---- */

  double* thomson = (double *)malloc(atmos.Nspace * sizeof(double));
  // Initialise the functions that use internal storage out here using the thomson
  // array as temporary storage
  cmo_Hminus_ff(500.0, thomson);
  cmo_H2minus_ff(500.0, thomson);
  cmo_H2plus_ff(500.0, thomson);
  Thomson(thomson);

  Atom* He = (atmos.elements[1].model) ? atmos.elements[1].model : NULL;
  

  for (int e = 0; e < atmos.Nelem; ++e)
  {
    Element* elem = &atmos.elements[e];
    if (!elem->n)
    {
      elem->n = matrix_double(elem->Nstage, atmos.Nspace);
      LTEpops_elem(elem);
    }
  }

  BackgroundState bg[input.Nthreads];
  for (int i = 0; i < input.Nthreads; ++i)
  {
    bg[i].chi_c = (double*)malloc(NrecStokes * atmos.Nspace * sizeof(double));
    bg[i].eta_c = (double*)malloc(NrecStokes * atmos.Nspace * sizeof(double));
    bg[i].sca_c = (double*)malloc(atmos.Nspace * sizeof(double));
    bg[i].Bnu = (double*)malloc(atmos.Nspace * sizeof(double));

    bg[i].chi = (double*)malloc(NrecStokes * atmos.Nspace * sizeof(double));
    bg[i].eta = (double*)malloc(NrecStokes * atmos.Nspace * sizeof(double));
    bg[i].scatt = (double*)malloc(atmos.Nspace * sizeof(double));

    if (atmos.Stokes && input.magneto_optical)
    {
      bg[i].chip = (double*)malloc(3 * atmos.Nspace * sizeof(double));
      bg[i].chip_c = (double*)malloc(3 * atmos.Nspace * sizeof(double));
    }
    else
    {
      bg[i].chip = NULL;
      bg[i].chip_c = NULL;
    }

    if (atmos.moving || atmos.Stokes)
    {
      bg[i].chi_ai = (double*)malloc(atmos.Nspace * sizeof(double));
      bg[i].eta_ai = (double*)malloc(atmos.Nspace * sizeof(double));
      bg[i].sca_ai = (double*)malloc(atmos.Nspace * sizeof(double));
    }
    else
    {
      bg[i].chi_ai = bg[i].chi_c;
      bg[i].eta_ai = bg[i].eta_c;
      bg[i].sca_ai = bg[i].sca_c;
    }

    bg[i].He = He;
    bg[i].NrecStokes = NrecStokes;
    bg[i].do_fudge = do_fudge;
    bg[i].Nfudge = Nfudge;
    bg[i].lambda_fudge = lambda_fudge;
    bg[i].fudge = fudge;
    bg[i].thomson = thomson;
    bg[i].bbState = (PassiveBbState*)malloc(sizeof(PassiveBbState));
    init_PassiveBbState(bg[i].bbState);
  }

  

  /* --- Allocate memory for the boolean array that stores whether
         a wavelength overlaps with a Bound-Bound transition in the
         background, or whether it is polarized --     -------------- */

  // atmos.backgrflags = (flags *)malloc(spectrum.Nspect * sizeof(flags));
  // for (int nspect = 0; nspect < spectrum.Nspect; nspect++) {
    // atmos.backgrflags[nspect].hasline = FALSE;
    // atmos.backgrflags[nspect].ispolarized = FALSE;
    // atmos.backgrflags[nspect] = 0;
  // }
  memset(atmos.backgrflags, 0, spectrum.Nspect * sizeof(flags));
  /* --- Allocate memory for the list of record numbers that specifies
         for each wavelength where to find the background opacity,
         scattering opacity, and emissivity --         -------------- */


  // if (atmos.moving || atmos.Stokes) {
  //   atmos.backgrrecno =
  //       (long *)calloc(2 * spectrum.Nspect * atmos.Nrays, sizeof(long));
  // } else
  //   atmos.backgrrecno = (long *)calloc(spectrum.Nspect, sizeof(long));
  if (atmos.moving || atmos.Stokes) {
        memset(atmos.backgrrecno, 0, 
               2 * spectrum.Nspect * atmos.Nrays * sizeof(long));
  } else
        memset(atmos.backgrrecno, 0, 
               spectrum.Nspect * sizeof(long));

  if (input.Nthreads > 1)
  {
    struct sched_task task;
    scheduler_add(&input.sched, &task, inner_background_sched, 
                  bg, spectrum.Nspect, input.workSize);
    scheduler_join(&input.sched, &task);
  }
  else
  {
    for (int nspect = 0; nspect < spectrum.Nspect; ++nspect)
    // for (int nspect = spectrum.Nspect - 1; nspect >= 0; --nspect)
    {
      inner_background(&bg[0], nspect, 0);
    }
  }
  

  if (atmos.moving || atmos.Stokes)
  {
    int maxIdx = (spectrum.Nspect - 1) * (atmos.Nrays * 2) + 
                 (atmos.Nrays - 1) * 2 + 2;
    int recAccumulate = 0;
    for (int i = 0; i < maxIdx; ++i)
    {
      if (atmos.backgrrecno[i] != 0)
      {
        atmos.backgrrecno[i] += recAccumulate;
        recAccumulate = atmos.backgrrecno[i];
      }
    }
    for (int i = maxIdx - 1; i >= 0; --i)
    {
      if (atmos.backgrrecno[i] == 0)
      {
        atmos.backgrrecno[i] = atmos.backgrrecno[i+1];
      }
    }
    for (int i = 1; i < maxIdx; ++i)
    {
      atmos.backgrrecno[i] -= atmos.backgrrecno[0];
      if (atmos.backgrrecno[i] < 0)
      {
        atmos.backgrrecno[i] = 0;
      }
    }
    atmos.backgrrecno[0] = 0;
  }
  else
  {
    for (int nspect = 0; nspect < spectrum.Nspect; ++nspect)
    {
      int lamu = nspect;
      if (lamu != 0)
      {
        atmos.backgrrecno[lamu] += atmos.backgrrecno[lamu-1];
      }
    }
    // NOTE(cmo): This shouldn't be necessary
    // for (int nspect = spectrum.Nspect - 1; nspect >= 0; --nspect)
    // {
    //   if (nspect != 0)
    //   {
    //     atmos.backgrrecno[nspect] = atmos.backgrrecno[nspect-1];
    //   }
    //   else
    //   {
    //     atmos.backgrrecno[nspect] = 0;
    //   }
    // }
  }

  // FILE* rec = fopen("bgrec.txt", "w");
  // char buf[512];
  // char comma[] = ",\n";
  // int maxIdx = spectrum.Nspect * (2 * atmos.Nrays) + atmos.Nrays*2 + 2;
  // for (int i = 0; i < maxIdx; ++i)
  // {
  //   snprintf(buf, 512, "%ld", atmos.backgrrecno[i]);
  //   fwrite(buf, 1, strlen(buf), rec);
  //   if (i != maxIdx - 1)
  //   {
  //     fwrite(comma, 1, strlen(comma), rec);
  //   }
  // }
  // fflush(rec);
  // fclose(rec);

  if (strcmp(input.background_File, "none"))
  {
    if ((atmos.fd_background =
            open(input.background_File, O_RDWR | O_CREAT, PERMISSIONS)) == -1) {
      sprintf(messageStr, "Unable to open output file %s", input.background_File);
      Error(ERROR_LEVEL_2, __func__, messageStr);
    }

    if (atmos.moving || atmos.Stokes)
    {
      long prevRec = -1;
      for (int nspect = 0; nspect < spectrum.Nspect; ++nspect)
      {
        for (int mu = 0; mu < atmos.Nrays; ++mu)
        {
          for (int to_obs = 0; to_obs <= 1; ++to_obs)
          {
            int i = nspect * (atmos.Nrays * 2) + mu * 2 + to_obs;
            if (atmos.backgrrecno[i] != prevRec)
            {
              // printf("%d, %d, %d, %d, %ld\n", nspect, mu, to_obs, i, atmos.backgrrecno[i]);
              writeBackground(nspect, mu, to_obs, atmos.chi_b[i], 
                              atmos.eta_b[i], atmos.sca_b[i], 
                              atmos.chip_b ? atmos.chip_b[i] : NULL);
              prevRec = atmos.backgrrecno[i];
            }
          }
        }
      }
    }
    else
    {
      for (int nspect = 0; nspect < spectrum.Nspect; ++nspect)
      {
        writeBackground(nspect, 0, 0, atmos.chi_b[nspect], 
                        atmos.eta_b[nspect], atmos.sca_b[nspect], 
                        NULL);
      }
    }


    if (write_analyze_output) {
      /* --- Write background record structure --          ------------ */

      writeBRS();

      /* --- Write out the metals and molecules --         ------------ */

      writeMetals("metals.out");
      writeMolecules(MOLECULAR_CONCENTRATION_FILE);
    }
  }

  
  // NOTE(cmo): Don't free background atoms, as we will need them again
  // if (atmos.Natom > 1) {
  //   for (int n = 1; n < atmos.Natom; n++)
  //     if (!atmos.atoms[n].active && !atmos.hydrostatic &&
  //         input.solve_ne < ITERATION)
  //       freeAtom(&atmos.atoms[n]);
  // }
  // if (atmos.Nmolecule > 1) {
  //   for (int n = 1; n < atmos.Nmolecule; n++)
  //     if (!atmos.molecules[n].active && !atmos.hydrostatic &&
  //         input.solve_ne < ITERATION)
  //       freeMolecule(&atmos.molecules[n]);
  // }

  // NOTE(cmo): Don't free the Kurucz data, we might need it again
  // if (strcmp(input.KuruczData, "none")) {
  //   free(atmos.Tpf);
  //   atmos.Tpf = NULL;
  //   for (int n = 0; n < atmos.Nelem; n++) {
  //     free(atmos.elements[n].ionpot);
  //     freeMatrix((void **)atmos.elements[n].pf);
  //     if (atmos.elements[n].n)
  //       freeMatrix((void **)atmos.elements[n].n);
  //   }
  // }

  cmo_Hminus_ff(0.0, NULL);
  cmo_H2minus_ff(0.0, NULL);
  cmo_H2plus_ff(0.0, NULL);

  free(thomson);
  for (int i = 0; i < input.Nthreads; ++i)
  {
    free(bg[i].chi);
    free(bg[i].eta);
    free(bg[i].scatt);

    free(bg[i].Bnu);

    free(bg[i].chi_c);
    free(bg[i].eta_c);
    free(bg[i].sca_c);

    free_PassiveBbState(bg[i].bbState);
    free(bg[i].bbState);

    if (atmos.moving || atmos.Stokes) {
      free(bg[i].chi_ai);
      free(bg[i].eta_ai);
      free(bg[i].sca_ai);
    }
    if (atmos.Stokes && input.magneto_optical) {
      free(bg[i].chip);
      free(bg[i].chip_c);
    }

    if (do_fudge) {
      free(lambda_fudge);
      freeMatrix((void **)fudge);
    }
  }
  // free(chi);
  // free(eta);
  // free(scatt);
  // free(Bnu);
  // free(thomson);
  // free(chi_c);
  // free(eta_c);
  // free(sca_c);

  // if (atmos.moving || atmos.Stokes) {
  //   free(chi_ai);
  //   free(eta_ai);
  //   free(sca_ai);
  // }
  // if (atmos.Stokes && input.magneto_optical) {
  //   free(chip);
  //   free(chip_c);
  // }

  // if (do_fudge) {
  //   free(lambda_fudge);
  //   freeMatrix((void **)fudge);
  // }

  CMO_PROF_FUNC_END();
}
/* ------- begin -------------------------- Background.c ------------ */

void Background(bool_t write_analyze_output, bool_t equilibria_only) {
  const char routineName[] = "Background";
  register int k, nspect, n, mu, to_obs;

  static int ne_iter = 0;
  char inputLine[MAX_LINE_SIZE];
  bool_t exit_on_EOF, do_fudge, fromscratch;
  int backgrrecno, index, Nfudge, NrecStokes;
  double *chi, *eta, *scatt, wavelength, *thomson, *chi_ai, *eta_ai, *sca_ai,
      Hmin_fudge, scatt_fudge, metal_fudge, *lambda_fudge, **fudge, *Bnu,
      *chi_c, *eta_c, *sca_c, *chip, *chip_c;
  Atom *He;
  FILE *fp_fudge;
  flags backgrflags;

  CMO_PROF_FUNC_START();

  getCPU(2, TIME_START, NULL);

  if (input.solve_ne == ONCE || input.solve_ne == ITERATION) {
    fromscratch = (input.solve_ne == ONCE ||
                   (input.solve_ne == ITERATION && ne_iter == 0))
                      ? TRUE
                      : FALSE;
    Solve_ne(atmos.ne, fromscratch);
    ne_iter++;
  }
  SetLTEQuantities();

  if (input.NonICE)
    readMolecules(MOLECULAR_CONCENTRATION_FILE);
  else
    ChemicalEquilibrium(N_MAX_CHEM_ITER, CHEM_ITER_LIMIT);

  if (equilibria_only) {

    /* --- If we only need ne, LTE populations and collisions, and
           chemical equilibrium leave here --          -------------- */

    getCPU(2, TIME_POLL, "Total Background");
    CMO_PROF_FUNC_END();
    return;
  }

  if (input.old_background) {
    if ((atmos.fd_background = open(input.background_File, O_RDONLY, 0)) ==
        -1) {
      sprintf(messageStr, "Unable to open input file %s",
              input.background_File);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    readBRS();
    CMO_PROF_FUNC_END();
    return;
  }

  getCPU(3, TIME_START, NULL);
  if (strcmp(input.fudgeData, "none")) {
    do_fudge = TRUE;

    /* --- Read wavelength-dependent fudge factors to compensate for
           missing UV backround line haze --           -------------- */

    if ((fp_fudge = fopen(input.fudgeData, "r")) == NULL) {
      sprintf(messageStr, "Unable to open input file %s", input.fudgeData);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    sprintf(messageStr, "\n-Fudging background opacities with file\n  %s\n\n",
            input.fudgeData);
    Error(MESSAGE, routineName, messageStr);

    getLine(fp_fudge, COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
    sscanf(inputLine, "%d", &Nfudge);
    lambda_fudge = (double *)malloc(Nfudge * sizeof(double));
    fudge = matrix_double(3, Nfudge);
    for (n = 0; n < Nfudge; n++) {
      getLine(fp_fudge, COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
      sscanf(inputLine, "%lf %lf %lf %lf", &lambda_fudge[n], &fudge[0][n],
             &fudge[1][n], &fudge[2][n]);
    }
    for (n = 0; n < 3 * Nfudge; n++)
      fudge[0][n] += 1.0;
    fclose(fp_fudge);
  } else
    do_fudge = FALSE;

  /* --- Allocate temporary storage space. The quantities are used
         for the following purposes:

       - chi, eta, scatt: Get contributions to opacity, emissivity,
         and scattering opacity, respectively, from a specific process
         for a given wavelength and possibly angle.

       - chi_c, eta_c, sca_c: Collect the total opacity, emissivity
         and scattering opacity for a given wavelength and possibly
         angle.

       - chi_ai, eta_ai: Collect the angle-independent part of
         opacity and emissivity for each wavelength so that these
         need not be recalculated in an angle-dependent case.
         When the atmosphere is not moving and has no magnetic fields
         these just point to the total quantities chi_c and eta_c.

   Note: In case of magnetic fields in the atmosphere chi, eta and
         chip, and chi_c, eta_c and chip_c contain all four Stokes
         parameters, and should be allocated a 4 and 3 times larger
         storage space, respectively.
         --                                            -------------- */

  if (atmos.Stokes)
    NrecStokes = 4;
  else
    NrecStokes = 1;

  chi_c = (double *)malloc(NrecStokes * atmos.Nspace * sizeof(double));
  eta_c = (double *)malloc(NrecStokes * atmos.Nspace * sizeof(double));
  sca_c = (double *)malloc(atmos.Nspace * sizeof(double));

  chi = (double *)malloc(NrecStokes * atmos.Nspace * sizeof(double));
  eta = (double *)malloc(NrecStokes * atmos.Nspace * sizeof(double));
  scatt = (double *)malloc(atmos.Nspace * sizeof(double));

  if (atmos.Stokes && input.magneto_optical) {
    chip = (double *)malloc(3 * atmos.Nspace * sizeof(double));
    chip_c = (double *)malloc(3 * atmos.Nspace * sizeof(double));
  } else {
    chip = NULL;
    chip_c = NULL;
  }

  if (atmos.moving || atmos.Stokes) {
    chi_ai = (double *)malloc(atmos.Nspace * sizeof(double));
    eta_ai = (double *)malloc(atmos.Nspace * sizeof(double));
    sca_ai = (double *)malloc(atmos.Nspace * sizeof(double));
  } else {
    chi_ai = chi_c;
    eta_ai = eta_c;
    sca_ai = sca_c;
  }
  Bnu = (double *)malloc(atmos.Nspace * sizeof(double));

  /* --- Thomson scattering by free electrons is wavelength independent
         in non-relativistic limit so we compute it only once -- ---- */

  thomson = (double *)malloc(atmos.Nspace * sizeof(double));
  Thomson(thomson);

  /* --- Check whether an atomic model is present for He -- --------- */

  He = (atmos.elements[1].model) ? atmos.elements[1].model : NULL;

  /* --- Read background files from Kurucz data file -- ------------- */

  atmos.Nrlk = 0;
  readKuruczLines(input.KuruczData);
  if (atmos.Nrlk > 0) {
    qsort(atmos.rlk_lines, atmos.Nrlk, sizeof(RLK_Line), rlk_ascend);
  }
  /* --- Allocate memory for the boolean array that stores whether
         a wavelength overlaps with a Bound-Bound transition in the
         background, or whether it is polarized --     -------------- */

  atmos.backgrflags = (flags *)malloc(spectrum.Nspect * sizeof(flags));
  for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
    // atmos.backgrflags[nspect].hasline = FALSE;
    // atmos.backgrflags[nspect].ispolarized = FALSE;
    atmos.backgrflags[nspect] = 0;
  }
  /* --- Allocate memory for the list of record numbers that specifies
         for each wavelength where to find the background opacity,
         scattering opacity, and emissivity --         -------------- */

  backgrrecno = 0;

  if (atmos.moving || atmos.Stokes) {
    atmos.backgrrecno =
        (long *)malloc(2 * spectrum.Nspect * atmos.Nrays * sizeof(long));
  } else
    atmos.backgrrecno = (long *)malloc(spectrum.Nspect * sizeof(long));

  /* --- Open output file for background opacity, emissivity,
         scattering --                                 -------------- */

  if ((atmos.fd_background =
           open(input.background_File, O_RDWR | O_CREAT, PERMISSIONS)) == -1) {
    sprintf(messageStr, "Unable to open output file %s", input.background_File);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  /* --- Go through the spectrum and add the different opacity and
         emissivity contributions. This is the main loop --  -------- */

  for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
    wavelength = spectrum.lambda[nspect];

    /* --- The Planck function at this wavelength --   -------------- */

    Planck(atmos.Nspace, atmos.T, wavelength, Bnu);

    /* --- Initialize the flags for this wavelength -- -------------- */

    // atmos.backgrflags[nspect].hasline = FALSE;
    // atmos.backgrflags[nspect].ispolarized = FALSE;
    atmos.backgrflags[nspect] = 0;

    /* --- Initialize angle-independent quantities --  -------------- */

    for (k = 0; k < atmos.Nspace; k++) {
      chi_ai[k] = 0.0;
      eta_ai[k] = 0.0;
      sca_ai[k] = thomson[k];
    }
    /* --- Negative hydrogen ion, bound-free and free-free -- ------- */

    if (Hminus_bf(wavelength, chi, eta)) {
      for (k = 0; k < atmos.Nspace; k++) {
        chi_ai[k] += chi[k];
        eta_ai[k] += eta[k];
      }
    }
    if (Hminus_ff(wavelength, chi)) {
      for (k = 0; k < atmos.Nspace; k++) {
        chi_ai[k] += chi[k];
        eta_ai[k] += chi[k] * Bnu[k];
      }
    }
    /* --- Opacity fudge factors, applied to Hminus opacity -- ------ */

    if (do_fudge) {
      Linear(Nfudge, lambda_fudge, fudge[0], 1, &wavelength, &Hmin_fudge,
             FALSE);
      for (k = 0; k < atmos.Nspace; k++) {
        chi_ai[k] *= Hmin_fudge;
        eta_ai[k] *= Hmin_fudge;
      }
    }
    /* --- Opacities from bound-free transitions in OH and CH -- ---- */

    if (OH_bf_opac(wavelength, chi, eta)) {
      for (k = 0; k < atmos.Nspace; k++) {
        chi_ai[k] += chi[k];
        eta_ai[k] += eta[k];
      }
    }
    if (CH_bf_opac(wavelength, chi, eta)) {
      for (k = 0; k < atmos.Nspace; k++) {
        chi_ai[k] += chi[k];
        eta_ai[k] += eta[k];
      }
    }
    /* --- Neutral Hydrogen Bound-Free and Free-Free --  ------------ */

    if (Hydrogen_bf(wavelength, chi, eta)) {
      for (k = 0; k < atmos.Nspace; k++) {
        chi_ai[k] += chi[k];
        eta_ai[k] += eta[k];
      }
    }
    Hydrogen_ff(wavelength, chi);
    for (k = 0; k < atmos.Nspace; k++) {
      chi_ai[k] += chi[k];
      eta_ai[k] += chi[k] * Bnu[k];
    }
    /* --- Rayleigh scattering by neutral hydrogen --  -------------- */

    if (Rayleigh(wavelength, atmos.H, scatt)) {
      for (k = 0; k < atmos.Nspace; k++) {
        sca_ai[k] += scatt[k];
      }
    }
    /* --- Rayleigh scattering by neutral helium --    -------------- */

    if (He && Rayleigh(wavelength, He, scatt)) {
      for (k = 0; k < atmos.Nspace; k++) {
        sca_ai[k] += scatt[k];
      }
    }
    /* --- Absorption by H + H^+ (referred to as H2plus free-free) -- */

    if (H2plus_ff(wavelength, chi)) {
      for (k = 0; k < atmos.Nspace; k++) {
        chi_ai[k] += chi[k];
        eta_ai[k] += chi[k] * Bnu[k];
      }
    }
    /* --- Rayleigh scattering and free-free absorption by
           molecular hydrogen --                       -------------- */

    if (Rayleigh_H2(wavelength, scatt)) {
      for (k = 0; k < atmos.Nspace; k++) {
        sca_ai[k] += scatt[k];
      }
    }
    if (H2minus_ff(wavelength, chi)) {
      for (k = 0; k < atmos.Nspace; k++) {
        chi_ai[k] += chi[k];
        eta_ai[k] += chi[k] * Bnu[k];
      }
    }
    /* --- Bound-Free opacities due to ``metals'' --   -------------- */

    if (do_fudge) {
      Linear(Nfudge, lambda_fudge, fudge[2], 1, &wavelength, &metal_fudge,
             FALSE);
    } else {
      metal_fudge = 1.0;
    }
    /* --- Note: Hydrogen bound-free opacities are calculated in
           routine Hydrogen_bf --                      -------------- */

    Metal_bf(wavelength, atmos.Natom - 1, atmos.atoms + 1, chi, eta);
    for (k = 0; k < atmos.Nspace; k++) {
      chi_ai[k] += chi[k] * metal_fudge;
      eta_ai[k] += eta[k] * metal_fudge;
    }
    /* --- Add the scattering opacity to the absorption part to store
           the total opacity --                        -------------- */

    if (do_fudge) {
      Linear(Nfudge, lambda_fudge, fudge[1], 1, &wavelength, &scatt_fudge,
             FALSE);
    } else {
      scatt_fudge = 1.0;
    }
    for (k = 0; k < atmos.Nspace; k++) {
      sca_ai[k] *= scatt_fudge;
      chi_ai[k] += sca_ai[k];
    }
    /* --- Now the contributions that may be angle-dependent due to the
           presence of atomic or molecular lines --    -------------- */

    if (atmos.moving || atmos.Stokes) {
      for (mu = 0; mu < atmos.Nrays; mu++) {
        for (to_obs = 0; to_obs <= 1; to_obs++) {
          index = 2 * (nspect * atmos.Nrays + mu) + to_obs;

          /* --- First, copy the angle-independent parts -- --------- */

          for (k = 0; k < atmos.Nspace; k++) {
            chi_c[k] = chi_ai[k];
            eta_c[k] = eta_ai[k];
            sca_c[k] = sca_ai[k];
          }

          /* --- Zero the polarized quantities, if necessary -- ----- */

          if (atmos.Stokes) {
            for (k = atmos.Nspace; k < 4 * atmos.Nspace; k++) {
              chi_c[k] = 0.0;
              eta_c[k] = 0.0;
            }
            if (input.magneto_optical)
              for (k = 0; k < 3 * atmos.Nspace; k++)
                chip_c[k] = 0.0;
          }
          /* --- Add opacity from passive atomic lines (including
                 hydrogen) --                          -------------- */

          if (input.allow_passive_bb) {
            backgrflags =
                passive_bb(wavelength, nspect, mu, to_obs, chi, eta, chip);
            // if (backgrflags.hasline) {
            if (backgrflags & HAS_LINE) {
              // atmos.backgrflags[nspect].hasline = TRUE;
              atmos.backgrflags[nspect] |= HAS_LINE;
              // if (backgrflags.ispolarized) {
              if (backgrflags & IS_POLARIZED) {
                NrecStokes = 4;
                // atmos.backgrflags[nspect].ispolarized = TRUE;
                atmos.backgrflags[nspect] |= IS_POLARIZED;
                if (input.magneto_optical) {
                  for (k = 0; k < 3 * atmos.Nspace; k++)
                    chip_c[k] += chip[k];
                }
              } else
                NrecStokes = 1;

              for (k = 0; k < NrecStokes * atmos.Nspace; k++) {
                chi_c[k] += chi[k];
                eta_c[k] += eta[k];
              }
            }
          }
          /* --- Add opacity from Kurucz line list --  -------------- */

          if (atmos.Nrlk > 0) {
            backgrflags = rlk_opacity(wavelength, nspect, mu, to_obs, chi, eta,
                                      scatt, chip);
            if (backgrflags & HAS_LINE) {
              // atmos.backgrflags[nspect].hasline = TRUE;
              atmos.backgrflags[nspect] |= HAS_LINE;
              // if (backgrflags.ispolarized) {
              if (backgrflags & IS_POLARIZED) {
                NrecStokes = 4;
                // atmos.backgrflags[nspect].ispolarized = TRUE;
                atmos.backgrflags[nspect] |= IS_POLARIZED;
                if (input.magneto_optical) {
                  for (k = 0; k < 3 * atmos.Nspace; k++)
                    chip_c[k] += chip[k];
                }
              } else
                NrecStokes = 1;

              for (k = 0; k < NrecStokes * atmos.Nspace; k++) {
                chi_c[k] += chi[k];
                eta_c[k] += eta[k];
              }
              if (input.rlkscatter) {
                for (k = 0; k < atmos.Nspace; k++) {
                  sca_c[k] += scatt[k];
                  chi_c[k] += scatt[k];
                }
              }
            }
          }
          /* --- Add opacity from molecular lines --   -------------- */

          backgrflags =
              MolecularOpacity(wavelength, nspect, mu, to_obs, chi, eta, chip);
          // if (backgrflags.hasline) {
          if (backgrflags & HAS_LINE) {
            // atmos.backgrflags[nspect].hasline = TRUE;
            atmos.backgrflags[nspect] |= HAS_LINE;
            // if (backgrflags.ispolarized) {
            if (backgrflags & IS_POLARIZED) {
              NrecStokes = 4;
              // atmos.backgrflags[nspect].ispolarized = TRUE;
              atmos.backgrflags[nspect] |= IS_POLARIZED;
              if (input.magneto_optical) {
                for (k = 0; k < 3 * atmos.Nspace; k++)
                  chip_c[k] += chip[k];
              }
            } else
              NrecStokes = 1;

            for (k = 0; k < NrecStokes * atmos.Nspace; k++) {
              chi_c[k] += chi[k];
              eta_c[k] += eta[k];
            }
          }
          /* --- Store angle-dependent results only if at least one
                 line was found at this wavelength --  -------------- */

          atmos.backgrrecno[index] = backgrrecno;
          if ((mu == atmos.Nrays - 1 && to_obs) ||
              // (atmos.backgrflags[nspect].hasline &&
              //  (atmos.moving || atmos.backgrflags[nspect].ispolarized))) {
              ((atmos.backgrflags[nspect] & HAS_LINE) &&
               (atmos.moving || (atmos.backgrflags[nspect] & IS_POLARIZED)))) {
            backgrrecno += writeBackground(nspect, mu, to_obs, chi_c, eta_c,
                                           sca_c, chip_c);
          }
        }
      }
    } else {
      /* --- Angle-independent case. First, add opacity from passive
             atomic lines (including hydrogen) --      -------------- */

      if (input.allow_passive_bb) {
        backgrflags = passive_bb(wavelength, nspect, 0, TRUE, chi, eta, NULL);
        // if (backgrflags.hasline) {
        if (backgrflags & HAS_LINE) {
          // atmos.backgrflags[nspect].hasline = TRUE;
          atmos.backgrflags[nspect] |= HAS_LINE;
          for (k = 0; k < atmos.Nspace; k++) {
            chi_c[k] += chi[k];
            eta_c[k] += eta[k];
          }
        }
      }
      /* --- Add opacity from Kurucz line list --      -------------- */

      if (atmos.Nrlk > 0) {
        backgrflags =
            rlk_opacity(wavelength, nspect, 0, TRUE, chi, eta, scatt, NULL);
        // if (backgrflags.hasline) {
        if (backgrflags & HAS_LINE) {
          // atmos.backgrflags[nspect].hasline = TRUE;
          atmos.backgrflags[nspect] |= HAS_LINE;
          for (k = 0; k < atmos.Nspace; k++) {
            chi_c[k] += chi[k];
            eta_c[k] += eta[k];
          }
          if (input.rlkscatter) {
            for (k = 0; k < atmos.Nspace; k++) {
              sca_c[k] += scatt[k];
              chi_c[k] += scatt[k];
            }
          }
        }
      }
      /* --- Add opacity from molecular lines --       -------------- */

      backgrflags =
          MolecularOpacity(wavelength, nspect, 0, TRUE, chi, eta, NULL);
      // if (backgrflags.hasline) {
      if (backgrflags & HAS_LINE) {
        // atmos.backgrflags[nspect].hasline = TRUE;
        atmos.backgrflags[nspect] |= HAS_LINE;
        for (k = 0; k < atmos.Nspace; k++) {
          chi_c[k] += chi[k];
          eta_c[k] += eta[k];
        }
      }
      /* --- Store results --                          -------------- */

      atmos.backgrrecno[nspect] = backgrrecno;
      backgrrecno += writeBackground(nspect, 0, 0, chi_c, eta_c, sca_c, NULL);
    }
  }

  if (write_analyze_output) {
    /* --- Write background record structure --          ------------ */

    writeBRS();

    /* --- Write out the metals and molecules --         ------------ */

    writeMetals("metals.out");
    writeMolecules(MOLECULAR_CONCENTRATION_FILE);
  }
  /* --- Clean up but keep H, H2, and active atom and/or molecule
         if appropriate --                               ------------ */

  if (atmos.Natom > 1) {
    for (n = 1; n < atmos.Natom; n++)
      if (!atmos.atoms[n].active && !atmos.hydrostatic &&
          input.solve_ne < ITERATION)
        freeAtom(&atmos.atoms[n]);
  }
  if (atmos.Nmolecule > 1) {
    for (n = 1; n < atmos.Nmolecule; n++)
      if (!atmos.molecules[n].active && !atmos.hydrostatic &&
          input.solve_ne < ITERATION)
        freeMolecule(&atmos.molecules[n]);
  }

  if (strcmp(input.KuruczData, "none")) {
    free(atmos.Tpf);
    atmos.Tpf = NULL;
    for (n = 0; n < atmos.Nelem; n++) {
      free(atmos.elements[n].ionpot);
      freeMatrix((void **)atmos.elements[n].pf);
      if (atmos.elements[n].n)
        freeMatrix((void **)atmos.elements[n].n);
    }
  }
  getCPU(3, TIME_POLL, "Background Opacity");

  // FILE* rec = fopen("bgrec_han.txt", "w");
  // char buf[512];
  // char comma[] = ",\n";
  // int maxIdx = spectrum.Nspect * (2 * atmos.Nrays) + atmos.Nrays*2 + 2;
  // for (int i = 0; i < maxIdx; ++i)
  // {
  //   snprintf(buf, 512, "%ld", atmos.backgrrecno[i]);
  //   fwrite(buf, 1, strlen(buf), rec);
  //   if (i != maxIdx - 1)
  //   {
  //     fwrite(comma, 1, strlen(comma), rec);
  //   }
  // }
  // fflush(rec);
  // fclose(rec);

  /* --- Free the temporary space allocated in the ff routines -- --- */

  Hminus_ff(0.0, NULL);
  H2minus_ff(0.0, NULL);
  H2plus_ff(0.0, NULL);

  free(chi);
  free(eta);
  free(scatt);
  free(Bnu);
  free(thomson);
  free(chi_c);
  free(eta_c);
  free(sca_c);

  if (atmos.moving || atmos.Stokes) {
    free(chi_ai);
    free(eta_ai);
    free(sca_ai);
  }
  if (atmos.Stokes && input.magneto_optical) {
    free(chip);
    free(chip_c);
  }

  if (do_fudge) {
    free(lambda_fudge);
    freeMatrix((void **)fudge);
  }
  getCPU(2, TIME_POLL, "Total Background");
  CMO_PROF_FUNC_END();
}
/* ------- end ---------------------------- Background.c ------------ */
