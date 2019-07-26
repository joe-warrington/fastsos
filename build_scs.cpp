// (c) 2019 ETH Zurich, Automatic Control Lab, Joe Warrington

#include <iostream>
#include <tuple>
#include <vector>
#include "polystring.h"

#include "timing.h"
#include "glbopts.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

#include "sos.h"

#include "build_scs.h"

using namespace std;

tuple<ScsCone *, ScsData *, ScsSolution *, ScsInfo> create_scs_model(
        PolyInfo& f_info, vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos,
        vector<vector<vector<int> > > &g_mono_exponents, vector<vector<vector<int> > > &h_mono_exponents,
        int dim, int deg, const vector<vector<int> > &exp_list_2d, int total_vars, vector<unsigned long int> &s_of_d_minus_djs,
        vector<unsigned long int> &s_of_d_minus_dj2s, string positivity_condition, int output_level) {

    timestamp_t t1, t2, t3, t4;

    if (output_level > 0) {
        cout << "\n  Creating model and determining number of variables and constraints..." << flush;
    }
    t1 = timenow();
    ScsCone *k = (ScsCone *) scs_calloc(1, sizeof(ScsCone));
    ScsData *d = (ScsData *) scs_calloc(1, sizeof(ScsData));
    ScsSolution *sol = (ScsSolution *) scs_calloc(1, sizeof(ScsSolution));
    ScsInfo info = {0};

    // Define problem dimensions
    int s_of_2d = (int) n_monomials(dim, 2 * deg);
    int s_of_d = (int) n_monomials(dim, deg);
    int psd_constr_rows = size_of_vmat(s_of_d);  // Number of elements in vectorized sigma_0 matrix
    vector<int> psd_side_lengths;
    vector<int> psd_constr_heights;
    psd_side_lengths.push_back(s_of_d);
    psd_constr_heights.push_back(psd_constr_rows);
    int n_g = g_infos.size();
    for (int j = 0; j < n_g; j++) {
        psd_side_lengths.push_back(s_of_d_minus_djs[j]);
        psd_constr_heights.push_back(size_of_vmat(s_of_d_minus_djs[j]));
        psd_constr_rows += psd_constr_heights[j + 1];
    }

    scs_int n = (scs_int) total_vars; // Length of optimization vector x
    scs_int m = (scs_int) (s_of_2d + psd_constr_rows); // Number of equality constraints in Ax + s = b
    t2 = timenow();
    if (output_level > 0) {
        cout << " done in " << time_string(t2 - t1) << ".\n" << flush;
        cout << "    SCS model has n = " << n << " variables (cols) and " << s_of_2d << " + " << psd_constr_rows
             << " = " << m << " constraints (rows).\n" << flush;
    }
    d->stgs = (ScsSettings *) scs_calloc(1, sizeof(ScsSettings));
    d->m = m;
    d->n = n;
    SCS(set_default_settings)(d);
    d->stgs->verbose = 0;
    d->stgs->acceleration_lookback = 10;
    d->stgs->rho_x = 1e-3;

    if (positivity_condition == "PSD") {
        // Create mixture of zero, free, and PSD cone variables
        k->f = (scs_int) s_of_2d;
//        cout << "  k->f = " << k->f << ", ";
        k->ssize = (scs_int) (n_g + 1);
//        cout << "k->ssize = " << k->ssize << ", k->s = {";
        scs_int *psd_side_lengths_array = new scs_int(psd_side_lengths.size());  // Allocate memory for s array
        copy(psd_side_lengths.begin(), psd_side_lengths.end(), psd_side_lengths_array);
        k->s = psd_side_lengths_array;
//        for (int j = 0; j < n_g + 1; j++) {
//            cout << k->s[j] << ", ";
//        }
//        cout << "\b\b}" << endl;
        k->l = 0;
        k->q = SCS_NULL;
        k->qsize = 0;
        k->ep = 0;
        k->ed = 0;
        k->p = SCS_NULL;
        k->psize = 0;
    } else {
        cout << "WARNING: PSD is the only positivity condition implemented with SCS." << endl;
    }

    ScsMatrix *A = d->A = (ScsMatrix *) scs_calloc(1, sizeof(ScsMatrix));
    scs_float *b = d->b = (scs_float *) scs_calloc(m, sizeof(scs_float));
    scs_float *c = d->c = (scs_float *) scs_calloc(n, sizeof(scs_float));

    // Compute number of nonzeros in A matrix before allocating memory to it
    if (output_level > 0) {
        cout << "  Computing number of nonzeros... " << flush;
    }
    t3 = timenow();
//    scs_int nnz = (scs_int) compute_scs_nonzeros_new(f_info, g_infos, h_infos, g_mono_exponents, h_mono_exponents,
//            dim, deg, exp_list_2d, s_of_d_minus_djs, s_of_d_minus_dj2s);
    scs_int nnz = (scs_int) compute_scs_nonzeros(f_info, g_infos, h_infos, g_mono_exponents, h_mono_exponents,
                                                     dim, deg, exp_list_2d);
    t4 = timenow();
    if (output_level > 0)
        cout << "Done in " << time_string(t4 - t3) << ". The A matrix (SCS standard form) has " << nnz << " nonzeros." << endl << flush;

    A->i = (scs_int *) scs_calloc(nnz, sizeof(scs_int));
    A->p = (scs_int *) scs_calloc((n + 1), sizeof(scs_int));
    A->x = (scs_float *) scs_calloc(nnz, sizeof(scs_float));
    A->n = d->n;
    A->m = d->m;

    return make_tuple(k, d, sol, info);
}

void create_scs_coeff_matches(tuple<ScsCone *, ScsData *, ScsSolution *, ScsInfo> M,
                                vector<double> &f_mono_coeffs, vector<vector<int> > &f_mono_exponents,
                                vector<vector<double> > &g_mono_coeffs, vector<vector<vector<int> > > &g_mono_exponents,
                                vector<vector<double> > &h_mono_coeffs, vector<vector<vector<int> > > &h_mono_exponents,
                                int dim, int deg, const vector<vector<int> > & exp_list_2d,
                                unsigned long int s_of_d, vector<unsigned long int> &s_of_d_minus_djs,
                                vector<unsigned long int> &s_of_d_minus_dj2s, int output_level) {

    ScsData *d = get<1>(M);

    ScsMatrix *A = d->A;
    scs_float *b = d->b;
    scs_float *c = d->c;

    int nnz_counter = 0;
    int col_counter = 0;
    int sigma_row_counter = 0;
    bool add_nz;
    // Update coefficients relating to lambda
    c[0] = -1.0;  // Corresponds to cost of lambda. "Minimize minus lambda"
    A->p[0] = 0; // Zeroth column starts at position 0
    A->i[0] = 0; // lambda appears in first column, in first position
    A->x[0] = 1.0; // lambda coefficient is 1 as we have "f - lambda = 0" when matching the coefficient of constants
    nnz_counter++; // update nnz_counter so that p[1] can be set correctly in the same manner as all later p[i]
    col_counter++;

    int s_of_2d = n_monomials(dim, deg * 2);
    unsigned long int idx;

    // Add entries for sigma_0
    int vmat_len = size_of_vmat(s_of_d);
    int k1 = -1, k2 = 0;
    vector<int> k1k2_exponent;
    vector<int> exp_to_match;
    for (int i = 0; i < vmat_len; i++) {  // Work left to right across A matrix
        A->p[i+col_counter] = nnz_counter; //cout << "A->p[" << i + col_counter << "] = " << nnz_counter << endl;
        k1++; // k1 starts at -1 before this loop, so that the first iteration yields k1 = 0.
        if (k1 == s_of_d) {
            k2++; // Increment row and column
            k1 = k2; // Set k1 to point to the diagonal entry of the matrix in column k2
        }
        // Add an entry for the coefficient matching
        k1k2_exponent = add_vecs(exp_list_2d[k1], exp_list_2d[k2]);
        auto it = find(exp_list_2d.begin(), exp_list_2d.end(), k1k2_exponent);
        if (it != exp_list_2d.end()) {
            idx = distance(exp_list_2d.begin(), it);
            A->i[nnz_counter] = idx; //cout << "A->i[" << nnz_counter << "] = " << A->i[nnz_counter] << endl;
            if (k1 == k2)
                A->x[nnz_counter] += 1.0;
            else
                A->x[nnz_counter] += sqrt(2);  // SCS stores the matrix with off-diagonals already scaled by sqrt(2)
            //cout << "A->x[" << nnz_counter << "] = " << A->x[nnz_counter] << endl;
            nnz_counter++;
        } else {cout << "WARNING: Coefficient match not found!\n";}
        // Add an entry to bind the element of sigma_0 to the PSD cone
        A->i[nnz_counter] = s_of_2d + i; //cout << "A->i[" << nnz_counter << "] = " << A->i[nnz_counter] << endl;
        A->x[nnz_counter] = -1.0; //cout << "A->x[" << nnz_counter << "] = " << A->x[nnz_counter] << endl;
        nnz_counter++;
    }
    col_counter += vmat_len;
    sigma_row_counter += vmat_len;

    // Loop through any g constraints present to populate constraints on corresponding sigma_j matrices
    int n_g = g_mono_coeffs.size();
    for (int j = 0; j < n_g; j++) {
        k1 = -1; k2 = 0;
        vmat_len = size_of_vmat(s_of_d_minus_djs[j]);
        for (int i = 0; i < vmat_len; i++) {  // Work left to right across A matrix
            A->p[i+col_counter] = nnz_counter; //cout << "A->p[" << i + col_counter << "] = " << nnz_counter << endl;
            k1++; // k1 starts at -1 before this loop, so that the first iteration yields k1 = 0.
            if (k1 == s_of_d_minus_djs[j]) {
                k2++; // Increment column
                k1 = k2; // Set k1 to point to the diagonal entry of the matrix in column k2
            }
            // Add an entry for the coefficient matching
            k1k2_exponent = add_vecs(exp_list_2d[k1], exp_list_2d[k2]);
            add_nz = false;
            for (int t = 0; t < g_mono_exponents[j].size(); t++) {
                exp_to_match = add_vecs(k1k2_exponent, g_mono_exponents[j][t]);
                auto it = find(exp_list_2d.begin(), exp_list_2d.end(), exp_to_match);
                if (it != exp_list_2d.end()) {
                    idx = distance(exp_list_2d.begin(), it);
                    A->i[nnz_counter] = idx;
                    //cout << "A->i[" << nnz_counter << "] = " << A->i[nnz_counter] << endl;
                    add_nz = A->x[nnz_counter] == 0;
                    if (k1 == k2)
                        A->x[nnz_counter] += g_mono_coeffs[j][t];
                    else
                        A->x[nnz_counter] += sqrt(2) * g_mono_coeffs[j][t];  // SCS stores the matrix with off-diagonals already scaled by sqrt(2)
                    //cout << "A->x[" << nnz_counter << "] = " << A->x[nnz_counter] << endl;
                } else { cout << "WARNING: Coefficient match not found!\n"; }
                if (add_nz) nnz_counter++;
            }
            // Add an entry to bind the element of sigma_j to the PSD cone
            A->i[nnz_counter] = s_of_2d + sigma_row_counter + i;
            //cout << "A->i[" << nnz_counter << "] = " << A->i[nnz_counter] << endl;
            A->x[nnz_counter] = -1.0;
            //cout << "A->x[" << nnz_counter << "] = " << A->x[nnz_counter] << endl;
            nnz_counter++;
        }
        col_counter += vmat_len;
        sigma_row_counter += vmat_len;
    }
    // Loop through any h constraints present to populate constraints on corresponding tau_j matrices
    int n_h = h_mono_coeffs.size();
    for (int j = 0; j < n_h; j++) {
        k1 = -1; k2 = 0;
        vmat_len = size_of_vmat(s_of_d_minus_dj2s[j]);
        for (int i = 0; i < vmat_len; i++) {  // Work left to right across A matrix
            A->p[i+col_counter] = nnz_counter; //cout << "A->p[" << i + col_counter << "] = " << nnz_counter << endl;
            k1++; // k1 starts at -1 before this loop, so that the first iteration yields k1 = 0.
            if (k1 == s_of_d_minus_dj2s[j]) {
                k2++; // Increment column
                k1 = k2; // Set k1 to point to the diagonal entry of the matrix in column k2
            }
            // Add an entry for the coefficient matching
            k1k2_exponent = add_vecs(exp_list_2d[k1], exp_list_2d[k2]);
            add_nz = false;
            for (int t = 0; t < h_mono_exponents[j].size(); t++) {
                exp_to_match = add_vecs(k1k2_exponent, h_mono_exponents[j][t]);
                auto it = find(exp_list_2d.begin(), exp_list_2d.end(), exp_to_match);
                if (it != exp_list_2d.end()) {
                    idx = distance(exp_list_2d.begin(), it);
                    A->i[nnz_counter] = idx;
                    //cout << "A->i[" << nnz_counter << "] = " << A->i[nnz_counter] << endl;
                    add_nz = A->x[nnz_counter] == 0;
                    if (k1 == k2)
                        A->x[nnz_counter] += h_mono_coeffs[j][t];
                    else
                        A->x[nnz_counter] += sqrt(2) * h_mono_coeffs[j][t];  // SCS stores the matrix with off-diagonals already scaled by sqrt(2)
                    //cout << "A->x[" << nnz_counter << "] = " << A->x[nnz_counter] << endl;
                } else { cout << "WARNING: Coefficient match not found!\n"; }
                if (add_nz) nnz_counter++;
            }
        }
        col_counter += vmat_len;
    }

    A->p[col_counter] = nnz_counter; //cout << "A->p[" << col_counter + 1 << "] = " << nnz_counter << endl;

    // Update RHS of equality constraint for coefficients of monomials in f

    for (int t = 0; t < f_mono_coeffs.size(); t++) {
        auto it = find(exp_list_2d.begin(), exp_list_2d.end(), f_mono_exponents[t]);
        if (it != exp_list_2d.end()) {
            idx = distance(exp_list_2d.begin(), it);
//            cout << "Found index in exponent list for the following term of objective f(x): ";
//            print_vec("", f_mono_exponents[t]);
            b[idx] += f_mono_coeffs[t];
        }
    }
}

int compute_scs_nonzeros_new(PolyInfo& f_info, vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos,
                             vector<vector<vector<int> > > & g_mono_exponents, vector<vector<vector<int> > > & h_mono_exponents,
                             int dim, int deg, const vector<vector<int> > & exp_list_2d, vector<unsigned long int> &s_of_d_minus_djs,
                             vector<unsigned long int> &s_of_d_minus_dj2s) {
    int nnz_counter = 0;
    int col_counter = 0;
    int sigma_row_counter = 0;
    bool add_nz;
    // Update count for lambda
    nnz_counter++; // update nnz_counter so that p[1] can be set correctly in the same manner as all later p[i]
    col_counter++;

    int s_of_d = n_monomials(dim, deg);
    int s_of_2d = n_monomials(dim, deg * 2);
    unsigned long int idx;

    // Add entries for sigma_0
    int vmat_len = size_of_vmat(s_of_d);
    int k1 = -1, k2 = 0;
    vector<int> k1k2_exponent;
    vector<int> exp_to_match;
    for (int i = 0; i < vmat_len; i++) {  // Work left to right across A matrix
        k1++; // k1 starts at -1 before this loop, so that the first iteration yields k1 = 0.
        if (k1 == s_of_d) {
            k2++; // Increment row and column
            k1 = k2; // Set k1 to point to the diagonal entry of the matrix in column k2
        }
        // Add an entry for the coefficient matching
        k1k2_exponent = add_vecs(exp_list_2d[k1], exp_list_2d[k2]);
        auto it = find(exp_list_2d.begin(), exp_list_2d.end(), k1k2_exponent);
        if (it != exp_list_2d.end()) {
            nnz_counter++;
        } else {cout << "WARNING: Coefficient match not found!\n";}
        // Add an entry to bind the element of sigma_0 to the PSD cone
        nnz_counter++;
    }
    col_counter += vmat_len;
    sigma_row_counter += vmat_len;

    // Loop through any g constraints present to populate constraints on corresponding sigma_j matrices
    int n_g = g_mono_exponents.size();
    for (int j = 0; j < n_g; j++) {
        k1 = -1; k2 = 0;
        vmat_len = size_of_vmat(s_of_d_minus_djs[j]);
        for (int i = 0; i < vmat_len; i++) {  // Work left to right across A matrix
            k1++; // k1 starts at -1 before this loop, so that the first iteration yields k1 = 0.
            if (k1 == s_of_d_minus_djs[j]) {
                k2++; // Increment column
                k1 = k2; // Set k1 to point to the diagonal entry of the matrix in column k2
            }
            // Add an entry for the coefficient matching
            k1k2_exponent = add_vecs(exp_list_2d[k1], exp_list_2d[k2]);
            add_nz = false;
            for (int t = 0; t < g_mono_exponents[j].size(); t++) {
                exp_to_match = add_vecs(k1k2_exponent, g_mono_exponents[j][t]);
                auto it = find(exp_list_2d.begin(), exp_list_2d.end(), exp_to_match);
                if (it != exp_list_2d.end()) {
                    add_nz = true;
                    break;
                } else { cout << "WARNING: Coefficient match not found!\n"; }
            }
            if (add_nz) nnz_counter++;
            // Add an entry to bind the element of sigma_j to the PSD cone
            nnz_counter++;
        }
        col_counter += vmat_len;
        sigma_row_counter += vmat_len;
    }
    // Loop through any h constraints present to populate constraints on corresponding tau_j matrices
    int n_h = h_mono_exponents.size();
    for (int j = 0; j < n_h; j++) {
        k1 = -1; k2 = 0;
        vmat_len = size_of_vmat(s_of_d_minus_dj2s[j]);
        for (int i = 0; i < vmat_len; i++) {  // Work left to right across A matrix
            k1++; // k1 starts at -1 before this loop, so that the first iteration yields k1 = 0.
            if (k1 == s_of_d_minus_dj2s[j]) {
                k2++; // Increment column
                k1 = k2; // Set k1 to point to the diagonal entry of the matrix in column k2
            }
            // Add an entry for the coefficient matching
            k1k2_exponent = add_vecs(exp_list_2d[k1], exp_list_2d[k2]);
            add_nz = false;
            for (int t = 0; t < h_mono_exponents[j].size(); t++) {
                exp_to_match = add_vecs(k1k2_exponent, h_mono_exponents[j][t]);
                auto it = find(exp_list_2d.begin(), exp_list_2d.end(), exp_to_match);
                if (it != exp_list_2d.end()) {
                    add_nz = true;
                    break;
                } else { cout << "WARNING: Coefficient match not found!\n"; }
            }
            if (add_nz) nnz_counter++;
        }
        col_counter += vmat_len;
    }
    return nnz_counter;
}

int compute_scs_nonzeros(PolyInfo& f_info, vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos,
        vector<vector<vector<int> > > & g_mono_exponents, vector<vector<vector<int> > > & h_mono_exponents,
        int dim, int deg, const vector<vector<int> > & exp_list_2d) {

    int nnz = 1;  // Accounts for the 1 in the (0, 0) entry of A, corresponding to lambda

    unsigned long int s_of_d = n_monomials(dim, deg);
    nnz += (s_of_d * (s_of_d + 1)) / 2; // Elements constraining sigma_0 to PSD cone.

    vector<unsigned long int> s_of_d_minus_djs;
    vector<unsigned long int> s_of_d_minus_dj2s;
    int dj;
    bool add_nz;
    for (int j = 0; j < g_infos.size(); j++) {
        // g_j(x) constraints cause nonzeros both in the upper (coefficient matching) and lower part (positivity) of A
        dj = ceil(g_infos[j].degree / 2.0);
        s_of_d_minus_djs.push_back(n_monomials(dim, deg - dj));
        nnz += (s_of_d_minus_djs[j] * (s_of_d_minus_djs[j] + 1)) / 2; // To allow for -1 entries for sigma_j in PSD
    }
    for (int j = 0; j < h_infos.size(); j++) {
        // h_j(x) constraints cause nonzeros only in the upper (coefficient matching) part of A
        dj = ceil(h_infos[j].degree / 2.0);
        s_of_d_minus_dj2s.push_back(n_monomials(dim, deg - dj));
    }
    for (unsigned long int i = 0; i < exp_list_2d.size(); i++) {
        // Match coefficients to see if there is a 1, sqrt(2), or 1/sqrt(2) appearing in the matrix.
        // This part loops through the rows of the A matrix, and then the columns

        // sigma_0 coefficient matching nonzeros
        for (int k2 = 0; k2 < s_of_d; k2++) {
            for (int k1 = k2; k1 < s_of_d; k1++) {
                if (exp_list_2d[i] == add_vecs(exp_list_2d[k1], exp_list_2d[k2]))
                    nnz += 1;
            }
        }
        // sigma_j coefficient matching nonzeros
        for (int j = 0; j < g_infos.size(); j++) {
            for (int k2 = 0; k2 < s_of_d_minus_djs[j]; k2++) {
                for (int k1 = k2; k1 < s_of_d_minus_djs[j]; k1++) {
                    add_nz = false;
                    for (int t = 0; t < g_infos[j].n_terms; t++) {
                        if (exp_list_2d[i] == add_vecs(g_mono_exponents[j][t],
                                                       add_vecs(exp_list_2d[k1], exp_list_2d[k2]))) {
                            add_nz = true;
                            break;
                        }
                    }
                    if (add_nz) nnz++;
                }
            }
        }
        // tau_j coefficient matching nonzeros
        for (int j = 0; j < h_infos.size(); j++) {
            for (int k2 = 0; k2 < s_of_d_minus_dj2s[j]; k2++) {
                for (int k1 = k2; k1 < s_of_d_minus_dj2s[j]; k1++) {
                    add_nz = false;
                    for (int t = 0; t < h_infos[j].n_terms; t++) {
                        if (exp_list_2d[i] == add_vecs(h_mono_exponents[j][t],
                                                       add_vecs(exp_list_2d[k1], exp_list_2d[k2]))) {
                            add_nz = true;
                            break;
                        }
                    }
                    if (add_nz) nnz++;
                }
            }
        }
    }
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

    int dj, side_length;

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

    // tau_j
    for (int j = 0; j < h_infos.size(); j++) {
        dj = ceil(h_infos[j].degree / 2.0);
        side_length = n_monomials(n, d - dj);
        mat_vsize = size_of_vmat(side_length);
        n_vars += mat_vsize;
        start_posns.push_back(vec_posn);
        labels.push_back("tau_" + to_string(j + 1));
        vec_posn += mat_vsize;
    }

    return make_tuple(n_vars, start_posns, labels);
}