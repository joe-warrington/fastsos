// (c) 2019 ETH Zurich, Automatic Control Lab, Joe Warrington

#include <iostream>
#include <iomanip>

#include "timing.h"
#include "manip_poly.h"
#include "build_mosek.h"
#include "build_scs.h"

#include <math.h>
#include <vector>
#include <tuple>

using namespace std;

const string COLOR_RED("\033[1;31m");
const string COLOR_YELLOW("\033[1;33m");
const string COLOR_GREEN("\033[1;32m");
const string COLOR_RESET("\033[0m");

unsigned long int a_choose_b(int a, int b){
    // Formula for a-choose-b is a!/[b!(a-b)!]
    if(a <= 0 || b <= 0) {
        cout << COLOR_RED << "Cannot compute a-choose-b for a = " << a << " and b = " << b << "!" << endl;
        exit(1);
    }
    else if (b > a) {
        cout << COLOR_RED << "Cannot compute a-choose-b for b greater than a! (" << b << " > " << a << ")" << endl;
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

tuple<double, int, string> sos_level_d(
        string& f_string, vector<string>& g_strings, vector<string>& h_strings,
        int d_request, string& positivity_condition, int output_level, string solver) {

    // 1. Parse polynomial strings and convert to tabular form
    auto data_tuple = create_polynomial_tables(f_string, g_strings, h_strings, output_level);
    // Update n and input data, removing unused dims
    int n = get<0>(data_tuple);
    vector<vector<int> > f_mono_exponents = get<7>(data_tuple);
    vector<vector<vector<int> > > g_mono_exponents = get<9>(data_tuple);
    vector<vector<vector<int> > > h_mono_exponents = get<11>(data_tuple);
    eliminate_unused_dims(n, f_mono_exponents, g_mono_exponents, h_mono_exponents, output_level);

    // 2. Work out minimum legal d for SOS problem and set d to this if necessary
    PolyInfo f_info = get<3>(data_tuple);
    vector<PolyInfo> g_infos = get<4>(data_tuple);
    vector<PolyInfo> h_infos = get<5>(data_tuple);
    int d = compute_legal_d(f_info, g_infos, h_infos, d_request);

    // 3. Pass to solver
    tuple<double, int, string> sol_tuple;
    if (solver == "mosek") {
        sol_tuple = solve_with_mosek(data_tuple, d, positivity_condition, output_level);
    } else if (solver == "scs") {
        sol_tuple = solve_with_scs(data_tuple, d, positivity_condition, output_level);
    } else {
        cout << "Unrecognized solver choice: " << solver << endl;
        sol_tuple = make_tuple(0, -99, "No solver");
    }

    return sol_tuple;
}