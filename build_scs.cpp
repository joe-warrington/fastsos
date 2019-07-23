//
// Created by Joe Warrington on 2019-07-23.
//

#include "build_scs.h"

#include <iostream>
#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

using namespace std;

static const char *simple_lp(void) {
    ScsCone *k = (ScsCone *) scs_calloc(1, sizeof(ScsCone));
    ScsData *d = (ScsData *) scs_calloc(1, sizeof(ScsData));
    ScsSolution *sol = (ScsSolution *) scs_calloc(1, sizeof(ScsSolution));
    ScsInfo info = {0};

    // Define problem dimensions
    scs_int n = 2; // Length of optimization vector x
    scs_int m = 3; // Number of equality constraints in Ax + s = b
    scs_int nnz = 4;  // Number of nonzeros in A matrix
    scs_int exitflag;
    scs_int success;

    d->stgs = (ScsSettings *) scs_calloc(1, sizeof(ScsSettings));
    d->m = m;
    d->n = n;
    SCS(set_default_settings)(d);

    k->f = 0;
    k->l = m - k->f;

    ScsMatrix *A = d->A = (ScsMatrix *) scs_calloc(1, sizeof(ScsMatrix));
    scs_float *b = d->b = (scs_float *) scs_calloc(m, sizeof(scs_float));
    scs_float *c = d->c = (scs_float *) scs_calloc(n, sizeof(scs_float));

    A->i = (scs_int *) scs_calloc(nnz, sizeof(scs_int));
    A->p = (scs_int *) scs_calloc((n + 1), sizeof(scs_int));
    A->x = (scs_float *) scs_calloc(nnz, sizeof(scs_float));
    A->n = d->n;
    A->m = d->m;

    // Populate cost function and constraint data
    c[0] = -1.0;
    c[1] = -1.0;
    b[0] = 0.0;
    b[1] = 0.0;
    b[2] = 2.4;

    A->p[0] = 0;
    A->p[1] = 2;
    A->p[2] = 4;
    A->i[0] = 0;
    A->i[1] = 2;
    A->i[2] = 1;
    A->i[3] = 2;
    A->x[0] = -1.0;
    A->x[1] = 1.0;
    A->x[2] = -1.0;
    A->x[3] = 1.0;

    /* this test fails with the default choice of 10 */
    d->stgs->acceleration_lookback = 10;

    bool ws = true;  // Generates segmentation fault if true...
    if (ws) {
        d->stgs->warm_start = 1;
        sol->x[0] = 1.2;
        sol->x[1] = 1.2;
        sol->s[0] = 1.2;
        sol->s[1] = 1.2;
        sol->s[2] = 0.0;
        sol->y[0] = 0.0;
        sol->y[1] = 0.0;
        sol->y[2] = 1.0;
    }

    // Solve LP
    exitflag = scs(d, k, sol, &info);

    // Collect results
    success = exitflag == SCS_SOLVED;

    cout << "x* = [";
    for (int i = 0; i < n; i++)
        cout << sol->x[i] << " ";
    cout << "\b]" << endl;
    cout << "y* = [";
    for (int i = 0; i < m; i++)
        cout << sol->y[i] << " ";
    cout << "\b]" << endl;
    cout << "s* = [";
    for (int i = 0; i < m; i++)
        cout << sol->s[i] << " ";
    cout << "\b]" << endl;

    SCS(free_data)(d, k);
    SCS(free_sol)(sol);
    mu_assert("joe_lp: SCS failed to produce outputflag SCS_SOLVED", success);
    return nullptr;
}