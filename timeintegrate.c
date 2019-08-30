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

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];

typedef struct
{
    double** nOld;
    double** nPrevIter;
} TimeDependentAtomicIterationState;

double
time_dependent_atom_update(double tStep, Atom* atom,
                           TimeDependentAtomicIterationState* atomState)
{
    int Nlevel = atom->Nlevel;

    double** A = matrix_double(Nlevel, Nlevel);
    double*  B = (double*)malloc(Nlevel * sizeof(double));
    double dPopsMax = 0.0;

    for (int k = 0; k < atmos.Nspace; ++k)
    {
        for (int i = 0; i < Nlevel; ++i)
        {
            B[i] = atomState->nOld[k][i];
            for (int j = 0; j < Nlevel; ++j)
            {
                int ij = i * Nlevel + j;
                int ji = j * Nlevel + i;
                A[i][j] = -atom->Gamma[ij][k] * tStep;
                A[j][i] = -atom->Gamma[ji][k] * tStep;
            }
        }
        for (int i = 0; i < Nlevel; ++i)
        {
            A[i][i] = 1.0;
            for (int j = 0; j < Nlevel; ++j)
            {
                if (j != i)
                {
                    A[i][i] -= A[j][i];
                }
            }
        }

        int solved = SolveLinearEq(Nlevel, A, B, TRUE);
        if (solved != 0)
            return -1.0;

        for (int i = 0; i < Nlevel; ++i)
        {
            atom->n[i][k] = B[i];
            double dPops = (B[i] / atomState->nPrevIter[k][i] - 1.0);
            atomState->nPrevIter[k][i] = B[i];

            dPopsMax = MAX(dPops, dPopsMax);
        }
    }
    // Need to add the Ng in here, rather than the time_dependent_iteration as would be done in the statEquil case.
    bool_t accel = Accelerate(atom->Ng_n, atom->n[0]);
    sprintf(messageStr, " %s,", atom->ID);
    double dPops = MaxChange(atom->Ng_n, messageStr, FALSE);
    Error(MESSAGE, NULL, (accel) ? " (accelerated)\n" : "\n");
    dPopsMax = MAX(dPops, dPopsMax);
    return dPopsMax;
}

int
time_dependent_pop_update(double tStep,
                          TimeDependentAtomicIterationState* atomState,
                          int nIter,
                          double iterLimit)
{
    double dPopsMax = 0.0;
    for (int nact = 0; nact < atmos.Nactiveatom; ++nact)
    {
        Atom* atom = atmos.activeatoms[nact];

        double dPops = time_dependent_atom_update(tStep, atom, 
                                                  &atomState[nact]);

        if (dPops < 0.0)
            return -1;
        // sprintf(messageStr, " %s, delta = %6.4E\n", atom->ID, dPops);
        // Error(MESSAGE, NULL, messageStr);

        dPopsMax = MAX(dPops, dPopsMax);
    }
    return dPopsMax < iterLimit;
}

int 
time_dependent_iteration(double tStep, int NmaxIter, double iterLimit)
{
    CMO_PROF_FUNC_START();
    bool_t evalOperator;

    if (NmaxIter <= 0)
    {
        CMO_PROF_FUNC_END();
        return 0;
    }

    double** JOld = matrix_double(spectrum.Nspect, atmos.Nspace);
    for (int n = 0; n < spectrum.Nspect; ++n)
    {
        for (int k = 0; k < atmos.Nspace; ++k)
        {
            double s = spectrum.J[n][k];
            if (isnan(s))
            {
                JOld[n][k] = 0.0;
            }
            else
            {
                JOld[n][k] = s;
            }
        }
    }

    for (int nact = 0; nact < atmos.Nactiveatom; ++nact)
    {
        Atom* atom = atmos.activeatoms[nact];
        atom->Ng_n = NgInit(atom->Nlevel * atmos.Nspace, input.Ngdelay,
                            input.Ngorder, input.Ngperiod, atom->n[0]);
    }

    TimeDependentAtomicIterationState activeAtomState[atmos.Nactiveatom];
    for (int nact = 0; nact < atmos.Nactiveatom; ++nact)
    {
        Atom* atom = atmos.activeatoms[nact];
        // nOld and nPrevIter are transposed relative to normal n, for speed of access
        // in the internal functions
        activeAtomState[nact].nOld = matrix_double(atmos.Nspace, atom->Nlevel);
        activeAtomState[nact].nPrevIter = matrix_double(atmos.Nspace, atom->Nlevel);
        for (int i = 0; i < atom->Nlevel; ++i)
        {
            for (int k = 0; k < atmos.Nspace; ++k)
            {
                activeAtomState[nact].nOld[k][i] 
                  = activeAtomState[nact].nPrevIter[k][i] 
                  = atom->n[i][k];
            }
        }
    }
    int niter = 1;

    while (niter <= NmaxIter)
    {
        for (int nact = 0; nact < atmos.Nactiveatom; ++nact)
        {
            initGammaAtom(atmos.activeatoms[nact], 1.0);
        }

        solve_spectrum_complete(evalOperator=TRUE);

        int converged = time_dependent_pop_update(tStep, activeAtomState, niter, iterLimit);

        if (converged == -1)
        {
            for (int nact = 0; nact < atmos.Nactiveatom; ++nact)
            {
                Atom* atom = atmos.activeatoms[nact];
                for (int i = 0; i < atom->Nlevel; ++i)
                {
                    for (int k = 0; k < atmos.Nspace; ++k)
                    {
                        atom->n[i][k] = activeAtomState[nact].nOld[k][i];
                    }
                }
            }
            for (int nspect = 0; nspect < spectrum.Nspect; ++nspect)
            {
                for (int k = 0; k < atmos.Nspace; ++k)
                {
                    spectrum.J[nspect][k] = JOld[nspect][k];
                }
            }
            for (int nact = 0; nact < atmos.Nactiveatom; ++nact)
            {
                Atom* a = atmos.activeatoms[nact];
                // freeMatrix((void**)a->Gamma);
                freeMatrix((void**)activeAtomState[nact].nOld);
                freeMatrix((void**)activeAtomState[nact].nPrevIter);
                NgFree(a->Ng_n);
            }
            freeMatrix((void**)JOld);
            return -1;
        }

        if (atmos.NPRDactive > 0)
        {
            double PRDiterlimit = input.PRDiterLimit;
            Redistribute(input.PRD_NmaxIter, PRDiterlimit);
        }

        if (converged)
            break;

        ++niter;

        if (input.solve_ne == ITERATION)
            cmo_background(FALSE, FALSE);
    }

    if (niter >= NmaxIter)
    {
        for (int nact = 0; nact < atmos.Nactiveatom; ++nact)
        {
            Atom* atom = atmos.activeatoms[nact];
            for (int i = 0; i < atom->Nlevel; ++i)
            {
                for (int k = 0; k < atmos.Nspace; ++k)
                {
                    atom->n[i][k] = activeAtomState[nact].nOld[k][i];
                }
            }
        }
        for (int nspect = 0; nspect < spectrum.Nspect; ++nspect)
        {
            for (int k = 0; k < atmos.Nspace; ++k)
            {
                spectrum.J[nspect][k] = JOld[nspect][k];
            }
        }
        for (int nact = 0; nact < atmos.Nactiveatom; ++nact)
        {
            Atom* a = atmos.activeatoms[nact];
            // freeMatrix((void**)a->Gamma);
            freeMatrix((void**)activeAtomState[nact].nOld);
            freeMatrix((void**)activeAtomState[nact].nPrevIter);
            NgFree(a->Ng_n);
        }
        freeMatrix((void**)JOld);
        return -1;
    }

    for (int nact = 0; nact < atmos.Nactiveatom; ++nact)
    {
        Atom* a = atmos.activeatoms[nact];
        // freeMatrix((void**)a->Gamma);
        freeMatrix((void**)activeAtomState[nact].nOld);
        freeMatrix((void**)activeAtomState[nact].nPrevIter);
        NgFree(a->Ng_n);
    }
    freeMatrix((void**)JOld);
    CMO_PROF_FUNC_END();
    return niter;
}