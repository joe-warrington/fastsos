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

tuple<ScsCone, ScsData, ScsSolution, ScsInfo> create_scs_model(PolyInfo& f_info, vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos,
                      int dim, int deg, int total_vars, string positivity_condition) {
    ScsCone *k = (ScsCone *) scs_calloc(1, sizeof(ScsCone));
    ScsData *d = (ScsData *) scs_calloc(1, sizeof(ScsData));
    ScsSolution *sol = (ScsSolution *) scs_calloc(1, sizeof(ScsSolution));
    ScsInfo info = {0};

    // Define problem dimensions
    int s_of_2d = (int) n_monomials(dim, 2 * deg);
    vector<int> psd_constr_heights(0, 0);
    int psd_constr_rows = 0;
    int n_g = g_infos.size();
    for (int j = 0; j < n_g; j++) {
        psd_constr_heights.push_back(size_of_vmat(n_monomials(dim, ceil(dim - g_infos[j].degree / 2.0))));
        psd_constr_rows += psd_constr_heights[j];
    }

    scs_int n = total_vars; // Length of optimization vector x
    scs_int m = s_of_2d + psd_constr_rows; // Number of equality constraints in Ax + s = b

    d->stgs = (ScsSettings *) scs_calloc(1, sizeof(ScsSettings));
    d->m = m;
    d->n = n;
    SCS(set_default_settings)(d);

    if (positivity_condition == "PSD") {
        // Create mixture of zero, free, and PSD cone variables
        k->f = s_of_2d;
        k->ssize = n_g;
        int psd_constr_heights_array[psd_constr_heights.size()];
        copy(psd_constr_heights.begin(), psd_constr_heights.end(), psd_constr_heights_array);
        k->s = psd_constr_heights_array;
    } else {
        cout << "PSD is the only positivity condition implemented with SCS." << endl;
    }

    ScsMatrix *A = d->A = (ScsMatrix *) scs_calloc(1, sizeof(ScsMatrix));
    scs_float *b = d->b = (scs_float *) scs_calloc(m, sizeof(scs_float));
    scs_float *c = d->c = (scs_float *) scs_calloc(n, sizeof(scs_float));

    // Compute number of nonzeros in A matrix before allocating memory to it
    vector<vector<int> > exp_list_2d = generate_all_exponents(dim, 2 * deg, 0);
    scs_int nnz = (scs_int) compute_scs_nonzeros(f_info, g_infos, h_infos, dim, deg, exp_list_2d);

    A->i = (scs_int *) scs_calloc(nnz, sizeof(scs_int));
    A->p = (scs_int *) scs_calloc((n + 1), sizeof(scs_int));
    A->x = (scs_float *) scs_calloc(nnz, sizeof(scs_float));
    A->n = d->n;
    A->m = d->m;

    // Populate cost function and constraint data
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


    return make_tuple(*k, *d, *sol, info);
}

void create_scs_coeff_matches(tuple<ScsCone, ScsData, ScsSolution> &M,
                                vector<double> &f_mono_coeffs, vector<vector<int> > &f_mono_exponents,
                                vector<vector<double> > &g_mono_coeffs, vector<vector<vector<int> > > &g_mono_exponents,
                                vector<vector<double> > &h_mono_coeffs, vector<vector<vector<int> > > &h_mono_exponents,
                                int dim, int deg, unsigned long int s_of_d, vector<unsigned long int> &s_of_d_minus_djs,
                                vector<unsigned long int> &s_of_d_minus_dj2s, int output_level) {

    ScsCone k = get<0>(M);
    ScsData d = get<1>(M);
    ScsSolution sol = get<2>(M);

    ScsMatrix *A = d.A;
    scs_float *b = d.b;
    scs_float *c = d.c;

    int nnz_counter = 0;
    // Update coefficients relating to lambda
    c[0] = -1.0;  // Corresponds to cost of lambda. "Minimize minus lambda"
    A->p[0] = 0; // Zeroth column starts at position 0
    A->i[0] = 0; // lambda appears in first column, in first position
    A->x[0] = 1.0; // lambda coefficient is 1 as we have "f - lambda = 0" when matching the coefficient of constants
    nnz_counter += 1; // update nnz_counter so that p[1] can be

    int s_of_2d = n_monomials(dim, deg * 2);
    vector<vector<int> > exp_list_2d = generate_all_exponents(dim, 2 * deg, 0);

    // Update RHS of equality constraint for coefficients of monomials in f
    unsigned long int idx;
    for (int t = 0; t < f_mono_coeffs.size(); t++) {
        auto it = find(exp_list_2d.begin(), exp_list_2d.end(), f_mono_exponents[t]);
        if (it != exp_list_2d.end()) {
            idx = distance(exp_list_2d.begin(), it);
            cout << "Found index in exponent list for the following term of objective f(x): ";
            print_vec("", f_mono_exponents[t]);
            b[idx] += f_mono_coeffs[t];
        }
    }


}

int compute_scs_nonzeros(PolyInfo& f_info, vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos,
        int dim, int deg, vector<vector<int> > exp_list_2d) {

    unsigned long int s_of_d = n_monomials(dim, deg);
    int nnz = (s_of_d * (s_of_d + 1)) / 2; // Elements constraining sigma_0 to PSD cone.

    vector<unsigned long int> s_of_d_minus_djs;
    vector<unsigned long int> s_of_d_minus_dj2s;
    int dj;
    for (int j = 0; j < g_infos.size(); j++) {
        dj = ceil(g_infos[j].degree / 2.0);
        s_of_d_minus_djs.push_back(n_monomials(dim, deg - dj));
        nnz += (s_of_d_minus_djs[j] * (s_of_d_minus_djs[j] + 1)) / 2; // To allow for -1 entries for sigma_j in PSD
    }
    for (int j = 0; j < h_infos.size(); j++) {
        dj = ceil(h_infos[j].degree / 2.0);
        s_of_d_minus_dj2s.push_back(n_monomials(dim, deg - dj));
    }
    for (unsigned long int i = 0; i < exp_list_2d.size(); i++) {
        // Match coefficients to see if there is a 1, sqrt(2), or 1/sqrt(2) appearing in the matrix.
        for (int k2 = 0; k2 < s_of_d; k2++) {
            for (int k1 = k2; k1 < s_of_d; k1++) {
                if (exp_list_2d[i] == add_vecs(exp_list_2d[k1], exp_list_2d[k2]))
                    nnz += 1;
            }
        }
    }
    cout << "The SCS A matrix has " << nnz << " nonzeros (not currently accounting for g and h eq constraints)." << endl;
    return nnz;
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