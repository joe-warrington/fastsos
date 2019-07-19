//
// Created by Joe Warrington on 2019-06-01.
//
#include <iostream>
#include <iomanip>
#include "mosek.h"
#include "fusion.h"
#include "polystring.h"
#include <sys/time.h>
#include <math.h>
#include <vector>
#include <tuple>

using namespace std;
using namespace mosek::fusion;
using namespace monty;

const string color_red("\033[1;31m");
const string color_yellow("\033[1;33m");
const string color_green("\033[1;32m");
const string color_reset("\033[0m");

typedef unsigned long long timestamp_t;
static timestamp_t timenow() {
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

string time_string(unsigned long long us_in) {
    ostringstream time_string_stream;     // Create an output string stream
    time_string_stream << std::fixed;    // Set Fixed -Point Notation
    time_string_stream << std::setprecision(3);    // Set precision to 2 digits

//    if (log10((double) us_in) < 3) {
//        time_string_stream << us_in;
//        return time_string_stream.str() + " us";
//    }
//    else
    if (log10((double) us_in) >= 6) {
        time_string_stream << us_in / 1000000.0;
        return time_string_stream.str() + " s";
    }
    else {
        time_string_stream << us_in / 1000.0;
        return time_string_stream.str() + " ms";
    }
}

unsigned long int a_choose_b(int a, int b){
    // Formula for a-choose-b is a!/[b!(a-b)!]
    if(a <= 0 || b <= 0) {
        cout << color_red << "Cannot compute a-choose-b for a = " << a << " and b = " << b << "!" << endl;
        exit(1);
    }
    else if (b > a) {
        cout << color_red << "Cannot compute a-choose-b for b greater than a! (" << b << " > " << a << ")" << endl;
        exit(1);
    }
    else if (b == a) {
        return 1;
    }
    else {  // Must have a > b > 0
        unsigned long int num = 1, den = 1;
        for (int i = b + 1; i <= a; i++) {
            num *= i;
        } // generates a! / b! for a > b, skipping all the cancelled terms 1 * 2 * ... * b
        for (int i = 1; i <= a - b; i++) {
            den *= i;
        } // generates (a - b)! in the denominator of the a-choose-b formula
        return num / den;
    }
}

unsigned long int n_monomials(int n, int d){
    unsigned long int result;
    int choosing_choice = max(n,d);
    // n+d choose d is the same as n+d choose n, but better overflow characteristics when the larger is chosen
    result = a_choose_b(n + d, choosing_choice);
    return result;
}

unsigned long int get_deg_start(int n, int deg_this_level) {
    // 0-indexed start index for monomials of degree deg_this_level for n-dimensional polynomials
    if (deg_this_level == 0) {return 0;}
    else if (deg_this_level == 1) {return 1;}
    else {
        // Start is just the number of monomials in a degree d-1 polynomial.
        return n_monomials(n, deg_this_level-1);
    }
}

template <class T>
void print_vec(const string& string_in, const T& v) {
    cout << string_in << ":\t";
    if (string_in.substr(0, 5) == "sigma") {  // Treat sigma (always matrices) differently, print in multiple rows
        cout << endl;
        int row_length = sqrt(v.size());
        for (int i = 0; i < row_length; i++) {
            for (int j = 0; j < row_length; j++) {
                cout << v[row_length * i + j] << "\t";
            }
            cout << endl;
        }
    }
    else {
        for (int i = 0; i < v.size(); i++)
            cout << v[i] << "\t";
        cout << endl;
    }
}

string monomial_to_string(const vector<int> v_in) {
    string s_out = "";
    for (int i = 0; i < v_in.size(); i++) {
        if (v_in[i] > 0)
            s_out += "x" + to_string(i + 1) + "^" + to_string(v_in[i]);
    }
    if (s_out == "")
        s_out = "1";  // If all exponents are zero, set the string representation to "1"
    return s_out;
}

vector<int> add_vecs(const vector <int> v1, const vector <int> v2) {
    vector<int> output(v1.size(), 0);
    for (int j = 0; j < v1.size(); j++)
        output[j] = v1[j] + v2[j];
    return output;
}

vector<int> sub_vecs(const vector <int> v1, const vector <int> v2) {
    vector<int> output(v1.size(), 0);
    for (int j = 0; j < v1.size(); j++)
        output[j] = v1[j] - v2[j];
    return output;
}

template <class T>
T sum_vec(const vector <T>& v) {
    T output = 0;
    for (int i = 0; i < v.size(); i++)
        output += v[i];
    return output;
}

void eliminate_unused_dims(int& n, vector<vector<int> >& f_exps, vector<vector<vector<int> > >& g_exps_list,
                           vector<vector<vector<int> > >& h_exps_list, int output_level) {
    // Go through dimensions and check that there is something non-zero in one of the constraints or the objective
    // Eliminate any dimensions where nothing appears.
    int m = g_exps_list.size();
    int p = h_exps_list.size();
    int dims_to_subtract = 0;
    bool dimension_used;

    for (int i = n - 1; i >= 0; i--) { // Loop over dimensions
        dimension_used = false;
        for (int t = 0; t < f_exps.size(); t++) {  // Loop over terms of f
            if (f_exps[t][i] != 0) {  // Check for non-zero entry in ith dimension of term t
                dimension_used = true;
                break;
            }
        }
        for (int j = 0; j < m; j++) {
            for (int t = 0; t < g_exps_list[j].size(); t++) {
                if (g_exps_list[j][t][i] != 0) {
                    dimension_used = true;
                    break;
                }
            }
        }
        for (int j = 0; j < p; j++) {
            for (int t = 0; t < h_exps_list[j].size(); t++) {
                if (h_exps_list[j][t][i] != 0) {
                    dimension_used = true;
                    break;
                }
            }
        }
        if (!dimension_used) {
            dims_to_subtract++;
            for (int t = 0; t < f_exps.size(); t++) {  // Loop over terms of f
                f_exps[t].erase(f_exps[t].begin() + i);
            }
            for (int j = 0; j < m; j++) {
                for (int t = 0; t < g_exps_list[j].size(); t++) {
                    g_exps_list[j][t].erase(g_exps_list[j][t].begin() + i);
                }
            }
            for (int j = 0; j < p; j++) {
                for (int t = 0; t < h_exps_list[j].size(); t++) {
                    h_exps_list[j][t].erase(h_exps_list[j][t].begin() + i);
                }
            }
            if (output_level > 0)
                cout << "Eliminated the unused dimension x" << i + 1 << " before any other processing took place.\n";
        }
    }

    n -= dims_to_subtract; // External scope will now see the reduced dimension of n, before calculating monomial lists

}

vector<vector <int> > generate_all_exponents(const int n, const int d, int output_level) {
    unsigned long int s_of_d = n_monomials(n, d);
    vector<int> blank_row(n, 0);
    auto row = blank_row;
    vector<vector <int> > vec_out(s_of_d, blank_row);
    timestamp_t t1, t2;
    if (output_level > 0)
        cout << "Generating all " << s_of_d << " exponents for n = " << n << " up to degree " << d << "...";
    t1 = timenow();
    int current_d = 1;
    vector<int> e_posns(d, 0);  // Positions of the up to d exponents
    bool new_d = true;  // Whether d was just incremented
    for (unsigned long int i = 1; i < s_of_d; i++) {
        // i = 0 corresponds to d = 0 and a row of all zeros, as initialized. Can therefore start from i = 1 (d = 1).
        row = blank_row;
        if (new_d) {
            for (int j = 0; j < current_d; j++)
                e_posns[j] = 0;  // Set position of all the current_d "ones" to 0.
                // Only the first (current_d) elements of e_posns are used.
            new_d = false;
        }
        else {
            if (e_posns[current_d - 1] == n - 1 && current_d > 1) {
                for (int j = current_d - 2; j >= 0; j--) {
                    // Search back through more significant "ones" to find one that hasn't reached the end posn (n - 1)
                    if (e_posns[j] < n - 1) {
                        e_posns[j] += 1;
                        for (int k = j + 1; k < current_d; k++)
                            e_posns[k] = e_posns[j];  // Set positions of all less significant "ones" to that of jth.
                        break;
                    }
                }
            } else {
                e_posns[current_d - 1] += 1;  // Advance least significant "one"
            }
        }
        for (int j = 0; j < current_d; j++)
            row[e_posns[j]] += 1;  // Add the exponents to the correct columns of the row
        vec_out[i] = row;  // Add the row to the output.

        if (row[n - 1] == current_d && i != s_of_d - 1) {
            new_d = true;
            current_d += 1;
        }
    }
    t2 = timenow();
    if (output_level > 0)
        cout << " done in " << time_string(t2 - t1) << "." << endl;
    return vec_out;
}

void constrain_to_cone(Model::t& M, Variable::t& matrix_var, const int matrix_size, const string& cone_type) {
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
    constrain_to_cone(M, sigma_0, s_of_d, positivity_condition);

    int m = g_infos.size();
    vector<Variable::t> sigma_j;
    vector<unsigned long int> s_of_d_minus_djs;
    int dj;
    for (int j = 0; j < m; j++) {
        dj = ceil(g_infos[j].degree / 2.0);
        s_of_d_minus_djs.push_back(n_monomials(n, d - dj));
        var_name = "sigma_" + to_string(j + 1);
        sigma_j.push_back(M->variable(var_name, Domain::unbounded((int) s_of_d_minus_djs[j], (int) s_of_d_minus_djs[j])));
        constrain_to_cone(M, sigma_j[j], s_of_d_minus_djs[j], positivity_condition);
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
        constrain_to_cone(M, tau_j[j], s_of_d_minus_dj2s[j], "Sym");
    }
    M->objective(ObjectiveSense::Maximize, lambda);
//    M->acceptedSolutionStatus(AccSolutionStatus.Anything);
//    M->setSolverParam("presolveUse", "off");
    return make_tuple(M, lambda, sigma_0, sigma_j, s_of_d_minus_djs, tau_j, s_of_d_minus_dj2s);
}

int compute_legal_d(PolyInfo f_info, vector<PolyInfo> g_infos, vector<PolyInfo> h_infos, int d_request) {
    int d_min = (int) ceil(f_info.degree / 2.0);
    for (int j = 0; j < g_infos.size(); j++)
        d_min = max(d_min, (int) ceil(g_infos[j].degree / 2.0));
    for (int j = 0; j < h_infos.size(); j++)
        d_min = max(d_min, (int) ceil(h_infos[j].degree / 2.0));
    if (d_request < d_min) {
        cout << "Using minimum legal d = " << d_min << ", which is > your choice of d = " << d_request << endl;
    }
    return max(d_min, d_request);
}

void create_coeff_matches(Model::t& M, vector<double>& f_mono_coeffs, vector<vector<int> >& f_mono_exponents,
                          vector<vector <double> >& g_mono_coeffs, vector<vector <vector <int> > >& g_mono_exponents,
                          vector<vector <double> >& h_mono_coeffs, vector<vector <vector <int> > >& h_mono_exponents,
                          int n, int d, unsigned long int s_of_d, Variable::t& lambda, Variable::t& sigma_0,
                          vector<Variable::t>& sigma_j, vector<unsigned long int>& s_of_d_minus_djs,
                          vector<Variable::t>& tau_j, vector<unsigned long int>& s_of_d_minus_dj2s, int output_level) {
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

tuple<double, ProblemStatus, SolutionStatus, SolutionStatus> sos_level_d(
        string& f_string, vector<string>& g_strings, vector<string>& h_strings,
        int d_request, string& positivity_condition, int output_level) {

    double obj_val = 0.0;
    ProblemStatus problem_status;
    SolutionStatus solution_status_primal;
    SolutionStatus solution_status_dual;

    // 0. State objective and constraints as strings read from file
    int m = g_strings.size();
    int p = h_strings.size();
    if (output_level > 0) {
        cout << "f(x)\t= " << f_string << endl;
        for (int i = 0; i < m; i++)
            cout << "g_" << i + 1 << "(x)\t= " << g_strings[i] << endl;
        for (int i = 0; i < p; i++)
            cout << "h_" << i + 1 << "(x)\t= " << h_strings[i] << endl;
    }
    // 1. Work out dimension and degree of each polynomial
    // Parse f(x)
    PolyInfo f_info = infer_dim_and_n_terms(f_string, false);
    vector<PolyInfo> g_infos;
    vector<PolyInfo> h_infos;
    for (int j = 0; j < m; j++)
        g_infos.push_back(infer_dim_and_n_terms(g_strings[j], false));
    for (int j = 0; j < p; j++)
        h_infos.push_back(infer_dim_and_n_terms(h_strings[j], false));
    int n = f_info.dimension;
    for (int j = 0; j < m; j++)
        n = max(n, g_infos[j].dimension);  // Increase n to match highest-dimensional ineq constraint if necessary
    for (int j = 0; j < p; j++)
        n = max(n, h_infos[j].dimension);  // Increase n to match highest-dimensional eq constraint if necessary
    f_info.dimension = n;
    for (int j = 0; j < m; j++)
        g_infos[j].dimension = n;
    for (int j = 0; j < p; j++)
        h_infos[j].dimension = n;

    if (f_info.degree == 0) {
        cout << "Cannot minimize a constant." << endl;
        return make_tuple(obj_val, problem_status, solution_status_primal, solution_status_dual);
    }
    vector<double> f_mono_coeffs(f_info.n_terms, 0.0); // List of monomial coefficients
    vector<vector<int> > f_mono_exponents(f_info.n_terms, vector<int>(f_info.dimension, 0));
    parse_poly(f_string, f_mono_coeffs, f_mono_exponents, f_info, false);

    // Parse g_1(x), ..., g_m(x)
    vector<vector <double> > g_mono_coeffs;
    vector<vector <vector <int> > > g_mono_exponents;
    for (int j = 0; j < m; j++) {
        if (g_infos.back().n_terms == 0) {
            g_infos.pop_back();  // Delete last polynomial because it was most likely just an empty space character
            cout << "  Polynomial g(" << j + 1 << ") has no terms: [" << g_strings[j] << "]. Skipping." << endl;
            continue;
        }
        vector<double> g_mono_coeffs_entry(g_infos[j].n_terms, 0.0);
        vector<vector<int> > g_mono_exponents_entry(g_infos[j].n_terms, vector<int>(g_infos[j].dimension, 0));
        g_mono_coeffs.push_back(g_mono_coeffs_entry);
        g_mono_exponents.push_back(g_mono_exponents_entry);
        parse_poly(g_strings[j], g_mono_coeffs[j], g_mono_exponents[j], g_infos[j], false);
    }
    m = g_infos.size(); // Update m to account for skipped strings containing no terms
    // Parse h_1(x), ..., h_p(x)
    vector<vector <double> > h_mono_coeffs;
    vector<vector <vector <int> > > h_mono_exponents;
    for (int j = 0; j < p; j++) {
        if (h_infos.back().n_terms == 0) {
            h_infos.pop_back();  // Delete last polynomial because it was most likely just an empty space character
            cout << "  Polynomial h(" << j + 1 << ") has no terms: [" << h_strings[j] << "]. Skipping." << endl;
            continue;
        }
        vector<double> h_mono_coeffs_entry(h_infos[j].n_terms, 0.0);
        vector<vector<int> > h_mono_exponents_entry(h_infos[j].n_terms, vector<int>(h_infos[j].dimension, 0));
        h_mono_coeffs.push_back(h_mono_coeffs_entry);
        h_mono_exponents.push_back(h_mono_exponents_entry);
        parse_poly(h_strings[j], h_mono_coeffs[j], h_mono_exponents[j], h_infos[j], false);
    }
    p = h_infos.size(); // Update m to account for skipped strings containing no terms

    eliminate_unused_dims(n, f_mono_exponents, g_mono_exponents, h_mono_exponents, 1);  // Update n and input data, removing unused dims

    // 2. Work out minimum legal d for SOS problem and set d to this if necessary
    int d = compute_legal_d(f_info, g_infos, h_infos, d_request);

    // 3. Work out number of PSD matrices and their sizes

    unsigned long int s_of_d = n_monomials(n, d);

    // 4. Create Mosek model and data entries of correct dimensions
    timestamp_t t1, t2;
    if (output_level > 0)
        cout << "Creating MOSEK variables and " << positivity_condition << " positivity constraints... " << flush;
    t1 = timenow();
        auto tuple_out = create_mosek_model(f_info, g_infos, h_infos, n, d, s_of_d, positivity_condition);
        Model::t M = get<0>(tuple_out); Variable::t lambda = get<1>(tuple_out);
        auto _M = finally([&]() { M->dispose(); });
        Variable::t sigma_0 = get<2>(tuple_out); vector<Variable::t> sigma_j = get<3>(tuple_out);
        vector<Variable::t> tau_j = get<5>(tuple_out);
        vector<unsigned long int> s_of_d_minus_djs = get<4>(tuple_out);
        vector<unsigned long int> s_of_d_minus_dj2s = get<6>(tuple_out);
    t2 = timenow();
    if (output_level > 0) {
        cout << "done in " << time_string(t2 - t1) << "." << endl << "  Matrix side lengths are " << s_of_d << ", ";
        for (int j = 0; j < m; j++) { cout << s_of_d_minus_djs[j] << ", "; }
        for (int j = 0; j < p; j++) { cout << s_of_d_minus_dj2s[j] << ", "; }
        cout << "\b\b." << endl;
    }
    // 5. Populate Mosek model from input data
    t1 = timenow();
        create_coeff_matches(M, f_mono_coeffs, f_mono_exponents, g_mono_coeffs, g_mono_exponents,
                             h_mono_coeffs, h_mono_exponents, n, d, s_of_d,
                             lambda, sigma_0, sigma_j, s_of_d_minus_djs, tau_j, s_of_d_minus_dj2s, output_level);
    t2 = timenow();
    if (output_level > 0)
        cout << "\b\b\b\b done in " << time_string(t2 - t1) << "." << endl;

    // 6. Solve and collect optimality information
    cout << "Solving..." << flush;
    try {
        t1 = timenow();
        M->solve();
        t2 = timenow();
        //Extract solution
        cout << " finished working in " << time_string(t2 - t1) << "." << endl;
        auto sol_lambda = lambda->level();
        cout << color_green << "  Lower bound for d = " << d << ": " << (*sol_lambda) << color_reset << endl;
        auto sol_sigma_0 = sigma_0->level();
        int max_matrix_print_size = 0;  // Change this hard-coded flag to print sigma matrices depending on size
        if (s_of_d <= max_matrix_print_size)
            print_vec("sigma_0", (*sol_sigma_0));
        for (int j = 0; j < m; j++) {
            auto sol_sigma_j = sigma_j[j]->level();
            if (s_of_d_minus_djs[j] <= max_matrix_print_size)
                print_vec("sigma_" + to_string(j + 1), (*sol_sigma_j));
        }
        for (int j = 0; j < p; j++) {
            auto sol_tau_j = tau_j[j]->level();
            if (s_of_d_minus_dj2s[j] <= max_matrix_print_size)
                print_vec("tau_" + to_string(j + 1), (*sol_tau_j));
        }
        obj_val = M->primalObjValue();
    }
    catch (mosek::fusion::SolutionError& e) {
        cout << color_yellow << " Didn't solve!\n  " << e.toString() << color_reset << endl;
    }
    catch (ParameterError& e) {
        cout << color_red << " Parameter error!\n  " << e.toString() << color_reset << endl;
    }
    catch (const exception& e) {
        cout << color_red << " Generic exception:\n  " << e.what() << color_reset << endl;
    }


    problem_status = M->getProblemStatus();
    solution_status_primal = M->getPrimalSolutionStatus();
    solution_status_dual = M->getDualSolutionStatus();

    return make_tuple(obj_val, problem_status, solution_status_primal, solution_status_dual);
}