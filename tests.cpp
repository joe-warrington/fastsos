//
// Created by Joe Warrington on 2019-06-27.
//

#include <iostream>
#include <vector>
#include "polystring.h"
#include "sos.h"

using namespace std;

const string color_red("\033[1;31m");
const string color_yellow("\033[1;33m");
const string color_green("\033[1;32m");
const string color_reset("\033[0m");

tuple<bool, bool, bool, PolyInfo, vector<double>, vector<vector<int> >> test_polynomial_parse(
        string s, vector<int> info_in, vector<double> coeffs_in, vector<vector<int> > exponents_in) {

    PolyInfo poly_info = infer_dim_and_n_terms(s, false);
    vector<double> coeffs_output(poly_info.n_terms, 0.0); // List of monomial coefficients
    vector<vector<int> > exponents_output(poly_info.n_terms, vector<int>(poly_info.dimension, 0));
    parse_poly(s, coeffs_output, exponents_output, poly_info, false);
    bool info_match = (info_in[0] == poly_info.dimension && info_in[1] == poly_info.degree && info_in[2] == poly_info.n_terms);
    bool coeffs_match = coeffs_in == coeffs_output;
    bool exponents_match = exponents_in == exponents_output;

    return make_tuple(info_match, coeffs_match, exponents_match, poly_info, coeffs_output, exponents_output);
}

void run_polynomial_parser_tests() {

    cout << "Testing polynomial parser..." << endl;

    vector<string> input_strings = {"1",
                                    "x1 + 1",
                                    "x1",
                                    "x1x2",
                                    "x1^1x2",
                                    "x1x2^1",
                                    "x1^1x2^1",
                                    "x1x2x1x1",
                                    "x1x2x1x1 + 2x1x2x1x2^1",
                                    "-3.5x1^1x2^1",
                                    "3.5x1^2x2^1",
                                    "3.5x1^2x3 - 2x1^2x3",
                                    "-3.5x1^2x3 + 2x1^2x3^1",
                                    "x1 + x3^2",
                                    "x1 + 1x3^2",
                                    "x1^1 + x2^2"};
    vector<vector<int> > info_to_match = {{1, 0, 1},  // Dimension, degree, number of terms
                                          {1, 1, 2},
                                          {1, 1, 1},
                                          {2, 2, 1},
                                          {2, 2, 1},
                                          {2, 2, 1},
                                          {2, 2, 1},
                                          {2, 4, 1},
                                          {2, 4, 2},
                                          {2, 2, 1},
                                          {2, 3, 1},
                                          {3, 3, 2},
                                          {3, 3, 2},
                                          {3, 2, 2},
                                          {3, 2, 2},
                                          {2, 2, 2}};
    vector<vector<double> > coeffs_to_match = {{1},  // Expected coefficients for each term
                                               {1, 1},
                                               {1},
                                               {1},
                                               {1},
                                               {1},
                                               {1},
                                               {1},
                                               {1, 2},
                                               {-3.5},
                                               {3.5},
                                               {3.5, -2},
                                               {-3.5, 2},
                                               {1, 1},
                                               {1, 1},
                                               {1, 1}};
    vector<vector<vector<int> > > exponents_to_match = {{{0}},  // Expected exponents for each term
                                                        {{1}, {0}},
                                                        {{1}},
                                                        {{1, 1}},
                                                        {{1, 1}},
                                                        {{1, 1}},
                                                        {{1, 1}},
                                                        {{3, 1}},
                                                        {{3, 1}, {2, 2}},
                                                        {{1, 1}},
                                                        {{2, 1}},
                                                        {{2, 0, 1}, {2, 0, 1}},
                                                        {{2, 0, 1}, {2, 0, 1}},
                                                        {{1, 0, 0}, {0, 0, 2}},
                                                        {{1, 0, 0}, {0, 0, 2}},
                                                        {{1, 0}, {0, 2}}};
    int poly_count = 0;
    int fail_count = 0;

    for (int i = 0; i < coeffs_to_match.size(); i++) {
        poly_count++;
        tuple<bool, bool, bool, PolyInfo, vector<double>, vector<vector<int> >> test_passed =
                test_polynomial_parse(input_strings[i], info_to_match[i], coeffs_to_match[i], exponents_to_match[i]);
        if (!get<0>(test_passed) || !get<1>(test_passed) || !get<2>(test_passed)) {
            cout << "#" << i + 1 << "\t" << color_red << "FAIL: ";
            if (!get<0>(test_passed)) // coeffs_match is false
                cout << "Polynomial info (dimension, degree, no. of terms) doesn't match for " << input_strings[i] << ". " << endl;
            if (!get<1>(test_passed)) // coeffs_match is false
                cout << "Coefficients don't match for " << input_strings[i] << ". " << endl;
            if (!get<2>(test_passed)) // exponents_match is false
                cout << "Exponents don't match for " << input_strings[i] << ". " << endl;
            if (!get<0>(test_passed)) {
                cout << "Expected polynomial info:" << endl;
                cout << "  Dimension " << info_to_match[i][0] << ", degree " << info_to_match[i][1] << ", n_terms " << info_to_match[i][2] << endl;
                cout << "Parsed polynomial info:" << endl;
                cout << "  Dimension " << get<3>(test_passed).dimension;
                cout << ", degree " << get<3>(test_passed).degree;
                cout << ", n_terms " << get<3>(test_passed).n_terms << endl;
            }
            if (!get<1>(test_passed) || !get<2>(test_passed)) {
                cout << "Expected result:" << endl;
                print_polynomial_table(coeffs_to_match[i], exponents_to_match[i]);
                cout << "Result:" << endl;
                print_polynomial_table(get<4>(test_passed), get<5>(test_passed));
            }
            fail_count++;
            cout << color_reset;
        } else
            cout  << "#" << i + 1 << "\t" << color_green << "PASS: " << input_strings[i] << " parsed correctly." << color_reset << endl;
    }
    if (fail_count > 0)
        cout << "WARNING: " << fail_count << "/" << poly_count << " polynomial parsing tests failed." << endl;
    else
        cout << "All " << poly_count << " polynomial parsing tests passed.\n";
        cout << "--------------------------------------------------------------\n\n";
}

void run_sos_program_tests() {
    cout << "Running SOS program tests..." << endl;
    vector<string> input_objs = {"x1",
                                 "x1",
                                 "x1",
                                 "x1 + x2"};
    vector<vector<string> > input_constrs = {{"x1", "1 - x1"},
                                             {},
                                             {"-x1", "x1 - 1"},
                                             {"-x_1^2 - x_2^2 + 1"}};
    vector<int> relaxation_degree = {2,
                                     2,
                                     2,
                                     2};
    vector<string> positivity_certs = {"PSD",
                                       "PSD",
                                       "PSD",
                                       "PSD"};
    vector<bool> primal_feasible = {true,
                                    true,
                                    false,
                                    true};
    vector<bool> primal_bounded = {true,
                                   false,
                                   true,
                                   true};
    vector<double> opt_values = {0.0,
                                 0.0,
                                 0.0,
                                 -1.41421};
    vector<double> sol_tols = {1e-6,
                               1e-6,
                               1e-6,
                               1e-4};

    int n_tests = input_objs.size();

    for (int i = 0; i < n_tests; i++) {
        sos_level_d(input_objs[i], input_constrs[i], relaxation_degree[i], positivity_certs[i]);
    }

}

void run_tests() {
    cout << "Running SOS code in test mode...\n";
    cout << "--------------------------------------------------------------\n\n";
    run_polynomial_parser_tests();
    run_sos_program_tests();
}
