//
// Created by Joe Warrington on 2019-07-23.
//


#include <iostream>
#include <tuple>
#include <vector>
#include "polystring.h"
#include "sos.h"
#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

#include "build_scs.h"


using namespace std;

tuple<ScsCone, ScsData, ScsSolution> create_scs_model(PolyInfo& f_info, vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos,
                      int dim, int deg, int total_vars, string positivity_condition) {
    ScsCone *k = (ScsCone *) scs_calloc(1, sizeof(ScsCone));
    ScsData *d = (ScsData *) scs_calloc(1, sizeof(ScsData));
    ScsSolution *sol = (ScsSolution *) scs_calloc(1, sizeof(ScsSolution));
    ScsInfo info = {0};

    // Define problem dimensions
    scs_int n = total_vars; // Length of optimization vector x
    scs_int m = 3; // Number of equality constraints in Ax + s = b
    scs_int nnz = 4;  // Number of nonzeros in A matrix
    scs_int exitflag;
    scs_int success;

    d->stgs = (ScsSettings *) scs_calloc(1, sizeof(ScsSettings));
    d->m = m;
    d->n = n;
    SCS(set_default_settings)(d);

    if (positivity_condition == "PSD") {
        // Create mixture of zero, free, and PSD cone variables
        k->f = 0;
        k->l = m - k->f;
    }

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

    return make_tuple(*k, *d, *sol);
}

int size_of_vmat(int side_length) {
    return (side_length * (side_length + 1)) / 2;
}

tuple<int, vector<int>, vector<string> > calc_n_vars(vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos,
        int n, int d) {
    int n_vars = 1; // Accounts for 'lambda'
    vector<int> start_posns(1, 0); // Start position of 0 for lambda
    vector<string> labels(1, "lambda");  // Label for lambda
    int vec_posn = 1; // Leave vector position at 1 to indicate start posn of next variable to be added

    // sigma_0
    unsigned long int s_of_d = n_monomials(n, d);
    int mat_vsize = size_of_vmat((int) s_of_d);
    n_vars += mat_vsize;
    start_posns.push_back(vec_posn);
    labels.push_back("sigma_0");
    vec_posn += mat_vsize;

    // tau_j
    int dj, side_length;
    for (int j = 0; j < h_infos.size(); j++) {
        dj = ceil(h_infos[j].degree / 2.0);
        side_length = n_monomials(n, d - dj);
        mat_vsize = size_of_vmat(side_length);
        n_vars += mat_vsize;
        start_posns.push_back(vec_posn);
        labels.push_back("tau_" + to_string(j + 1));
        vec_posn += mat_vsize;
    }

    // sigma_j
    for (int j = 0; j < g_infos.size(); j++) {
        dj = ceil(g_infos[j].degree / 2.0);
        side_length = n_monomials(n, d - dj);
        mat_vsize = size_of_vmat(side_length);
        n_vars += mat_vsize;
        start_posns.push_back(vec_posn);
        labels.push_back("sigma_" + to_string(j + 1));
        vec_posn += mat_vsize;
    }

    return make_tuple(n_vars, start_posns, labels);
}