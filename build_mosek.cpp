//
// Created by Joe Warrington on 2019-07-23.
//

#include <iostream>
#include <iomanip>
#include "mosek.h"
#include "fusion.h"
#include "polystring.h"
#include "build_mosek.h"
#include "sos.h"
#include <sys/time.h>
#include <math.h>
#include <vector>
#include <tuple>

using namespace std;
using namespace mosek::fusion;
using namespace monty;

void mosek_constrain_to_cone(Model::t &M, Variable::t &matrix_var, const int matrix_size, const string &cone_type) {
    if (cone_type == "Sym") {
        M->constraint(Expr::sub(matrix_var, Expr::transpose(matrix_var)), Domain::equalsTo(0.0));  // make symmetric
        return;
    } else if (cone_type == "PSD") {
        M->constraint(Expr::sub(matrix_var, Expr::transpose(matrix_var)), Domain::equalsTo(0.0));  // make symmetric
        M->constraint(matrix_var, Domain::inPSDCone(matrix_size));  // Symmetric positive definite matrix
        return;
    }
    else if (cone_type == "DSOS") {
        M->constraint(Expr::sub(matrix_var, Expr::transpose(matrix_var)), Domain::equalsTo(0.0));  // make symmetric
        Variable::t matrix_var_t = M->variable(Domain::unbounded((int) matrix_size, (int) matrix_size));
        M->constraint(Expr::sub(matrix_var_t, matrix_var), Domain::greaterThan(0.0));  // t >= matrix_var
        M->constraint(Expr::sub(matrix_var_t, Expr::mul(-1.0, matrix_var)), Domain::greaterThan(0.0)); // t >= -matrix_var
        M->constraint(Expr::mul(1.0, matrix_var->diag()), Domain::greaterThan(0.0));
        if (matrix_size > 1) {
            for (unsigned long int i = 0; i < matrix_size; i++) {
                auto cexpr = Expr::mul(1.0, matrix_var->index(new_array_ptr<int, 1>({(int) i, (int) i})));
                for (unsigned long int j = 0; j < matrix_size; j++) {
                    if (i != j)
                        cexpr = Expr::sub(cexpr, matrix_var_t->index(new_array_ptr<int, 1>({(int) i, (int) j})));
                }
                M->constraint(cexpr, Domain::greaterThan(0.0));
            }
        }
        return;
    } else if (cone_type == "SDSOS") {
        if (matrix_size <= 2) {
            M->constraint(matrix_var, Domain::inPSDCone(matrix_size));
        }
        else {
            M->constraint(Expr::sub(matrix_var, Expr::transpose(matrix_var)), Domain::equalsTo(0.0));  // make symmetric
            vector<Variable::t> small_mats;
            vector<Variable::t> small_mats_cone_elems;
            vector<Expression::t> blank_expr_row(matrix_size, Expr::constTerm(0.0));
            vector<vector<Expression::t> > constr_exprs(matrix_size, blank_expr_row);
            unsigned long int m_idx = 0;
            for (unsigned long int i = 0; i < matrix_size - 1; i++) {
                for (unsigned long int j = i + 1; j < matrix_size; j++) {
                    try {
                        small_mats.push_back(M->variable(Domain::unbounded(3)));
                        small_mats_cone_elems.push_back(M->variable(Domain::inQCone(3)));
                        M->constraint(Expr::sub(small_mats_cone_elems[m_idx]->index(0),
                                                Expr::add(small_mats[m_idx]->index(0), small_mats[m_idx]->index(2))), Domain::equalsTo(0.0));
                        M->constraint(Expr::sub(small_mats_cone_elems[m_idx]->index(1),
                                                Expr::mul(2.0, small_mats[m_idx]->index(1))), Domain::equalsTo(0.0));
                        M->constraint(Expr::sub(small_mats_cone_elems[m_idx]->index(2),
                                                Expr::sub(small_mats[m_idx]->index(0), small_mats[m_idx]->index(2))), Domain::equalsTo(0.0));
                        constr_exprs[i][i] = Expr::add(constr_exprs[i][i], small_mats[m_idx]->index(0));
                        constr_exprs[i][j] = Expr::add(constr_exprs[i][j], small_mats[m_idx]->index(1));
                        constr_exprs[j][j] = Expr::add(constr_exprs[j][j], small_mats[m_idx]->index(2));
//                        constr_exprs[i][i] = Expr::add(constr_exprs[i][i], Expr::mul(0.5, Expr::add(small_mats_cone_elems[m_idx]->index(0), small_mats_cone_elems[m_idx]->index(2))));
//                        constr_exprs[i][j] = Expr::add(constr_exprs[i][j], Expr::mul(0.5, small_mats_cone_elems[m_idx]->index(1)));
//                        constr_exprs[j][j] = Expr::add(constr_exprs[j][j], Expr::mul(0.5, Expr::sub(small_mats_cone_elems[m_idx]->index(0), small_mats_cone_elems[m_idx]->index(2))));
                        m_idx++;
                    }
                    catch (mosek::fusion::IndexError& e) {
                        cout << "i = " << i << ", j = " << j << ", matrix_size = " << matrix_size << endl;
                        cout << e.what() << endl;
                    }
                }
            }
            for (unsigned long int i = 0; i < matrix_size; i++) {
                for (unsigned long int j = i; j < matrix_size; j++) {
//                    cout << i << ", " << j << ": " << constr_exprs[i][j]->toString() << endl;
                    M->constraint(Expr::sub(matrix_var->index(i, j), constr_exprs[i][j]), Domain::equalsTo(0.0));
                    // Each matrix entry is equal to the sum of contributions from 2x2 submatrices
                }
            }
        }
        return;
    }
    else {
        cout << "Unknown positivity condition requested: " << cone_type << "." << endl;
        throw exception();
    }
}

tuple<Model::t, Variable::t, Variable::t, vector<Variable::t>, vector<unsigned long int>,
vector<Variable::t>, vector<unsigned long int> > create_mosek_model(
        PolyInfo& f_info, vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos,
        int n, int d, unsigned long int s_of_d, string positivity_condition) {
    Model::t M = new Model("SOS");
    // Create variables and objective function
    string var_name;
    Variable::t lambda = M->variable("lambda", Domain::unbounded());  // Scalar
    Variable::t sigma_0 = M->variable("sigma_0", Domain::unbounded((int) s_of_d, (int) s_of_d));  // Define dense matrix
    mosek_constrain_to_cone(M, sigma_0, s_of_d, positivity_condition);

    int m = g_infos.size();
    vector<Variable::t> sigma_j;
    vector<unsigned long int> s_of_d_minus_djs;
    int dj;
    for (int j = 0; j < m; j++) {
        dj = ceil(g_infos[j].degree / 2.0);
        s_of_d_minus_djs.push_back(n_monomials(n, d - dj));
        var_name = "sigma_" + to_string(j + 1);
        sigma_j.push_back(M->variable(var_name, Domain::unbounded((int) s_of_d_minus_djs[j], (int) s_of_d_minus_djs[j])));
        mosek_constrain_to_cone(M, sigma_j[j], s_of_d_minus_djs[j], positivity_condition);
    }
    int p = h_infos.size();
    vector<Variable::t> tau_j;
    vector<unsigned long int> s_of_d_minus_dj2s;
    int dj2;
    for (int j = 0; j < p; j++) {
        dj2 = ceil(h_infos[j].degree / 2.0);
        s_of_d_minus_dj2s.push_back(n_monomials(n, d - dj2));
        var_name = "tau_" + to_string(j + 1);
        tau_j.push_back(M->variable(var_name, Domain::unbounded((int) s_of_d_minus_dj2s[j], (int) s_of_d_minus_dj2s[j])));
        mosek_constrain_to_cone(M, tau_j[j], s_of_d_minus_dj2s[j], "Sym");
    }
    M->objective(ObjectiveSense::Maximize, lambda);
//    M->acceptedSolutionStatus(AccSolutionStatus.Anything);
//    M->setSolverParam("presolveUse", "off");
    return make_tuple(M, lambda, sigma_0, sigma_j, s_of_d_minus_djs, tau_j, s_of_d_minus_dj2s);
}


void create_mosek_coeff_matches(Model::t &M, vector<double> &f_mono_coeffs, vector<vector<int> > &f_mono_exponents,
                                vector<vector<double> > &g_mono_coeffs, vector<vector<vector<int> > > &g_mono_exponents,
vector<vector<double> > &h_mono_coeffs, vector<vector<vector<int> > > &h_mono_exponents,
int n, int d, unsigned long int s_of_d, Variable::t &lambda, Variable::t &sigma_0,
        vector<Variable::t> &sigma_j, vector<unsigned long int> &s_of_d_minus_djs,
        vector<Variable::t> &tau_j, vector<unsigned long int> &s_of_d_minus_dj2s,
int output_level) {
// For monomials in s(2d), find all entries (alpha, beta) for which x^(alpha+beta) = monomial.
// For the entries of s(2d) for which there is a term in f, set RHS of constraint to the coefficient. Otherwise 0.

unsigned long int s_of_2d = n_monomials(n, 2 * d);
vector<unsigned long int> exponent_degree_start_rows(2 * d + 2, 0);
for (int i = 1; i < 2 * d + 2; i++)
exponent_degree_start_rows[i] = n_monomials(n, i - 1);

int n_contributors;
double constr_rhs;
vector<vector <int> > all_2d_exponents = generate_all_exponents(n, 2 * d, output_level);
vector<vector <int> > all_d_exponents = all_2d_exponents; // Technically wasteful as it duplicates the 2d list
all_d_exponents.resize(s_of_d);  // Only retain the exponents up to degree d
int new_mod = 0, old_mod = 0;
auto constr_lhs = Expr::constTerm(0.0);
vector<int> exponent_to_add(n, 0);
//    int degree_to_add;
bool check_this_k1;

//    cout << "Eliminating unused monomials by looking for zero entries on the matrix diagonal... " << endl;
//    vector<vector<int> > reduced_d_exponents = all_d_exponents;
//    bool keep_this_monomial, present_in_this_term;
//    unsigned long int filtered_s_of_d = s_of_d;
//    for (unsigned long int i = s_of_d - 1; i >= 0; i--) {
//        keep_this_monomial = false;
//        // Iterate over terms in f
//        for (int t = 0; t < f_mono_coeffs.size(); t++) {
//            present_in_this_term = true;
//            for (int dim = 0; dim < n; d++) {
//                if (2 * all_d_exponents[i][dim] != f_mono_exponents[t][dim]) {
//                    present_in_this_term = false;
//                    break;
//                }
//            }
//            if (present_in_this_term) {
//                keep_this_monomial = true;
//                break;
//            }
//
//        }
//        if (keep_this_monomial) {
//            continue;
//        }
//        for (int j = 0; j < g_mono_coeffs.size(); j++) {  // for all constraint functions g_j(x)
//            for (int t = 0; t < g_mono_coeffs[j].size(); t++) {  // for all terms t in g_j(x)
//                present_in_this_term = true;
//                for (int k = 0; k < s_of_d_minus_djs[j]; k++) {  // for all multiplying terms in sigma_j(x)
//                    for (int dim = 0; dim < n; d++) {
//                        if (2 * all_d_exponents[i][dim] != g_mono_exponents[j][t][dim] + all_d_exponents[k][dim]) {
//                            present_in_this_term = false;
//                            break;
//                        }
//                    }
//                }
//                if (present_in_this_term) {
//                    keep_this_monomial = true;
//                    break; // break out of term loop
//                }
//            }
//            if (present_in_this_term) {
//                keep_this_monomial = true;
//                break;  // break out of j loop
//            }
//        }
//        if (!keep_this_monomial) {
//            cout << "  Deleted monomial  " << monomial_to_string(all_d_exponents[i]) << endl;
//            reduced_d_exponents.erase(i);
//        }
//    }
//    filtered_s_of_d = all_d_exponents.size();
if (output_level > 0) {
cout << "Creating the " << s_of_2d << " coefficient matching constraints... ";
cout << flush;
}
for (unsigned long int i = 0; i < s_of_2d; i++){
// Print progress in percent based on fraction of i indices covered.
old_mod = new_mod;
new_mod = (int) (100 * i) / s_of_2d;
if (new_mod > old_mod && output_level > 0) {
if (new_mod > 1)
cout << "\b\b\b\b" << setw(2) << setfill(' ') << new_mod << "% " << flush;
else
cout << setw(2) << setfill(' ') << new_mod << "% " << flush;
}
// For each monomial represented in the sigma_0 matrix, find row, col addresses of all terms contributing
n_contributors = 0;
if (i == 0)
constr_lhs = Expr::mul(1.0, lambda);  // Include lambda in the constant-term coefficient matching
else
constr_lhs = Expr::constTerm(0.0);
// Entries of sigma_0 matrix
for (unsigned long int k1 = 0; k1 < s_of_d; k1++) {
exponent_to_add = sub_vecs(all_2d_exponents[i], all_d_exponents[k1]);
check_this_k1 = true;
for (int l = 0; l < n; l++) {if (exponent_to_add[l] < 0) check_this_k1 = false;}
if (check_this_k1) {
//                degree_to_add = sum_vec(exponent_to_add);
//                for (unsigned long int k2 = max(k1, exponent_degree_start_rows[degree_to_add]);
//                     k2 < exponent_degree_start_rows[degree_to_add + 1]; k2++) {
for (unsigned long int k2 = k1; k2 < s_of_d; k2++) {
if (all_d_exponents[k2] == exponent_to_add) {
n_contributors++;
if (k1 == k2)
constr_lhs = Expr::add(constr_lhs,
                       Expr::mul(1.0, sigma_0->index(new_array_ptr<int, 1>({(int) k1, (int) k2}))));
else
constr_lhs = Expr::add(constr_lhs,
                       Expr::mul(2.0, sigma_0->index(new_array_ptr<int, 1>({(int) k1, (int) k2}))));
break;
}
}
}
}
for (int j = 0; j < g_mono_coeffs.size(); j++) {  // for all constraint functions g_j(x)
for (int t = 0; t < g_mono_coeffs[j].size(); t++) {  // for all terms t in g_j(x)
for (unsigned long int k1 = 0; k1 < s_of_d_minus_djs[j]; k1++) {
exponent_to_add = sub_vecs(all_2d_exponents[i],
                           add_vecs(all_d_exponents[k1], g_mono_exponents[j][t]));
check_this_k1 = true;
for (int l = 0; l < n; l++) {if (exponent_to_add[l] < 0) check_this_k1 = false;}
if (check_this_k1) {
//                        degree_to_add = sum_vec(exponent_to_add);
//                        for (unsigned long int k2 = max(k1, exponent_degree_start_rows[degree_to_add]);
//                             k2 < min(s_of_d_minus_djs[j], exponent_degree_start_rows[degree_to_add + 1]); k2++) {
for (unsigned long int k2 = k1; k2 < s_of_d_minus_djs[j]; k2++) {
if (all_d_exponents[k2] == exponent_to_add) {
if (k1 == k2)
constr_lhs = Expr::add(constr_lhs,
                       Expr::mul(1.0 * g_mono_coeffs[j][t],
                                 sigma_j[j]->index(new_array_ptr<int, 1>({(int) k1, (int) k2}))));
else
constr_lhs = Expr::add(constr_lhs,
                       Expr::mul(2.0 * g_mono_coeffs[j][t],
                                 sigma_j[j]->index(new_array_ptr<int, 1>({(int) k1, (int) k2}))));
break;
}
}
}
}
}
}
for (int j = 0; j < h_mono_coeffs.size(); j++) {  // for all constraint functions g_j(x)
for (int t = 0; t < h_mono_coeffs[j].size(); t++) {  // for all terms t in g_j(x)
for (unsigned long int k1 = 0; k1 < s_of_d_minus_dj2s[j]; k1++) {
exponent_to_add = sub_vecs(all_2d_exponents[i],
                           add_vecs(all_d_exponents[k1], h_mono_exponents[j][t]));
check_this_k1 = true;
for (int l = 0; l < n; l++) {if (exponent_to_add[l] < 0) check_this_k1 = false;}
if (check_this_k1) {
//                        degree_to_add = sum_vec(exponent_to_add);
//                        for (unsigned long int k2 = max(k1, exponent_degree_start_rows[degree_to_add]);
//                             k2 < min(s_of_d_minus_dj2s[j], exponent_degree_start_rows[degree_to_add + 1]); k2++) {
for (unsigned long int k2 = k1; k2 < s_of_d_minus_dj2s[j]; k2++) {
if (all_d_exponents[k2] == exponent_to_add) {
if (k1 == k2)
constr_lhs = Expr::add(constr_lhs,
                       Expr::mul(1.0 * h_mono_coeffs[j][t],
                                 tau_j[j]->index(new_array_ptr<int, 1>({(int) k1, (int) k2}))));
else
constr_lhs = Expr::add(constr_lhs,
                       Expr::mul(2.0 * h_mono_coeffs[j][t],
                                 tau_j[j]->index(new_array_ptr<int, 1>({(int) k1, (int) k2}))));
break;
}
}
}
}
}
}

constr_rhs = 0.0;
for (int t = 0; t < f_mono_coeffs.size(); t++) {
if (f_mono_exponents[t] == all_2d_exponents[i])
constr_rhs += f_mono_coeffs[t];
}
M->constraint(constr_lhs, Domain::equalsTo(constr_rhs));
}
}
