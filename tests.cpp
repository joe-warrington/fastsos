// (c) 2019 ETH Zurich, Automatic Control Lab, Joe Warrington

#include <iostream>
#include <vector>
#include "manip_poly.h"
#include "sos.h"

using namespace std;

const string COLOR_RED("\033[1;31m");
const string COLOR_YELLOW("\033[1;33m");
const string COLOR_GREEN("\033[1;32m");
const string COLOR_RESET("\033[0m");

tuple<bool, bool, bool, PolyInfo, vector<double>, vector<vector<int> >> test_polynomial_parse(
        string s, vector<int> info_in, vector<double> coeffs_in, vector<vector<int> > exponents_in) {

    PolyInfo poly_info = infer_dim_and_n_terms(s, false);
    vector<double> coeffs_output(poly_info.n_terms, 0.0); // List of monomial coefficients
    vector<vector<int> > exponents_output(poly_info.n_terms, vector<int>(poly_info.dimension, 0));
    parse_poly(s, coeffs_output, exponents_output, poly_info, false);
    bool info_match = (info_in[0] == poly_info.dimension &&
                       info_in[1] == poly_info.degree &&
                       info_in[2] == poly_info.n_terms);
    bool coeffs_match = coeffs_in == coeffs_output;
    bool exponents_match = exponents_in == exponents_output;

    return make_tuple(info_match, coeffs_match, exponents_match, poly_info, coeffs_output, exponents_output);
}

tuple<int, int> run_polynomial_parser_tests() {

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
    int test_count = 0;
    int fail_count = 0;

    for (int i = 0; i < coeffs_to_match.size(); i++) {
        test_count++;
        tuple<bool, bool, bool, PolyInfo, vector<double>, vector<vector<int> >> test_passed =
                test_polynomial_parse(input_strings[i], info_to_match[i], coeffs_to_match[i], exponents_to_match[i]);
        if (!get<0>(test_passed) || !get<1>(test_passed) || !get<2>(test_passed)) {
            cout << "#" << i + 1 << "\t" << COLOR_RED << "FAIL: ";
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
            cout << COLOR_RESET;
        } else
            cout  << "#" << i + 1 << "\t" << COLOR_GREEN << "PASS: " << input_strings[i] << " parsed correctly." << COLOR_RESET << endl;
    }
    if (fail_count > 0)
        cout << "WARNING: " << fail_count << "/" << test_count << " polynomial parsing tests failed." << endl;
    else
        cout << "All " << test_count << " polynomial parsing tests passed.\n";

    cout << "--------------------------------------------------------------\n\n";
    return make_tuple(test_count - fail_count, test_count);
}

tuple<bool, bool, double, int> test_sos_program(
        string& obj_in, vector<string> ineq_constrs_in, vector<string> eq_constrs_in, int d_in, string pos_cert_in,
        tuple<int, double, double, bool> expected_result_in, string solver_choice) {


    // Solve SOS problem

    int testing_output_level = 0;  // Controls level of text output during tests

    auto sol_info = sos_level_d(obj_in, ineq_constrs_in, eq_constrs_in, d_in, pos_cert_in,
                                testing_output_level, solver_choice);

    // Collect result and check against expected behaviour
    double obj_val = get<0>(sol_info);
    int sol_status = get<1>(sol_info);
    string sol_status_string = get<2>(sol_info);

    int expected_sol_status = get<0>(expected_result_in);
    double expected_obj_val = get<1>(expected_result_in);
    double expected_tolerance = get<2>(expected_result_in);
    bool expected_num_problems = get<3>(expected_result_in);

    bool obj_val_correct = false;
    bool report_obj_value = false;

    bool sol_status_correct = (sol_status == expected_sol_status);

    if (sol_status == 1 || sol_status == 2) {
        report_obj_value = true;
        if (abs(obj_val - expected_obj_val) <= expected_tolerance) {
            obj_val_correct = true;
        }
    }

    if (expected_sol_status != 1) {
        obj_val_correct = true;
        report_obj_value = false;
    }

    if (expected_num_problems) {
        obj_val_correct = true;
        report_obj_value = false;
        sol_status_correct = true; // Basically let the test pass if we know it causes numerical problems
    }

    if (report_obj_value) cout << "Objective value:\t" << obj_val << endl;
    cout << "Solution status: \t" << sol_status << " - " << sol_status_string << endl;

    return make_tuple(obj_val_correct, sol_status_correct, obj_val, sol_status);
}

tuple<int, int> run_sos_program_tests(string solver_choice) {
    cout << "Running SOS program tests..." << endl;
    vector<string> input_objs = {"x1^2 + 1.5",
                                 "x1",
                                 "x1",
                                 "0.5 + x1^1",
                                 "x2",
                                 "x1 + x2",
                                 "x1 + x2"};
    vector<vector<string> > input_ineq_constrs = {{},
                                                  {"-x1"},
                                                  {"-x1", "x1 - 1"},
                                                  {"-x1^2 - x2^2 + 1"},
                                                  {"-x2^2 - x5^2 + 1"},
                                                  {"-x1^2 + 7.0x1 - x2^2 + 7.0x2 - 12.25"},
                                                  {"-x1^2 + 7.0x1 - x2^2 + 7.0x2 - 12.25"}};
    vector<vector<string> > input_eq_constrs = {{},
                                                {},
                                                {},
                                                {},
                                                {},
                                                {"x1 - x3 - 2x4 - 4x5", "x2 - x6 - 2x7 - 4x8", "x3^2 - x3",
                                                 "x4^2 - x4", "x5^2 - x5", "x6^2 - x6", "x7^2 - x7", "x8^2 - x8"},
                                                {"x1 - x3 - 2x4 - 4x5", "x2 - x6 - 2x7 - 4x8", "x3^2 - x3",
                                                 "x4^2 - x4", "x5^2 - x5", "x6^2 - x6", "x7^2 - x7", "x8^2 - x8"}};
    vector<int> relaxation_degree = {1,
                                     2,
                                     2,
                                     2,
                                     2,
                                     1,
                                     3};
    vector<string> positivity_certs = {"PSD",
                                       "PSD",
                                       "PSD",
                                       "PSD",
                                       "PSD",
                                       "PSD",
                                       "PSD"};
    vector<int> sol_statuses = {1, // optimal solution found
                                -1, // unbounded
                                -2, // infeasible
                                1,
                                1,
                                1,
                                1};
    vector<double> opt_values = {1.5,
                                 0.0,
                                 0.0,
                                 -0.5,
                                 -1.0,
                                 2.05025,
                                 2.33647};
    vector<double> sol_tols = {1e-6,
                               1e-6,
                               1e-6,
                               1e-4,
                               1e-4,
                               1e-4,
                               1e-4};
    vector<bool> num_problems = {false,
                                 false,
                                 false,
                                 false,
                                 false,
                                 false,
                                 false};

    int n_tests = input_objs.size();

    int test_count = 0;
    int fail_count = 0;
    tuple<int, double, double, bool> expected_result;
    tuple<bool, bool, double, int> test_output; // Objective value correct, solver status correct, objective value

    for (int i = 0; i < n_tests; i++) {
        test_count++;
        string problem_string = "min " + input_objs[i];
        if (!input_ineq_constrs[i].empty() || !input_eq_constrs[i].empty()) problem_string += "   s.t.  ";
        for (int j = 0; j < input_ineq_constrs[i].size(); j++)
            problem_string += input_ineq_constrs[i][j] + " >= 0,  ";
        for (int j = 0; j < input_eq_constrs[i].size(); j++)
            problem_string += input_eq_constrs[i][j] + " = 0, ";
        cout << "#" << i + 1 << "\t" << problem_string << "\b\b " << endl;
        expected_result = make_tuple(sol_statuses[i], opt_values[i], sol_tols[i], num_problems[i]);

        test_output = test_sos_program(input_objs[i], input_ineq_constrs[i], input_eq_constrs[i],
                                       relaxation_degree[i], positivity_certs[i], expected_result, solver_choice);
        //test_output = (obj_val_correct, sol_status_correct, obj_val, sol_status)
        if (!get<0>(test_output) || !get<1>(test_output)) { // either objective value or solution status incorrect
            cout << COLOR_RED << "FAIL: ";
            if (!get<0>(test_output) && (get<3>(test_output) == 1 || get<3>(test_output) == 2)) {
                cout << "Objective value is incorrect: " << get<2>(test_output) << " instead of "
                        << opt_values[i] << " +- " << sol_tols[i] << endl;
            }
            if (!get<1>(test_output)) {
                cout << "Solution status is incorrect: " << get<3>(test_output) << " instead of " << sol_statuses[i] << endl;
            }
            fail_count++;
            cout << COLOR_RESET << "\n";
        } else
            cout << COLOR_GREEN << "PASS: solution behaved as expected.\n" << COLOR_RESET << endl;
    }
    if (fail_count > 0)
        cout << "WARNING: " << fail_count << "/" << test_count << " SOS programming tests failed." << endl;
    else
        cout << "All " << test_count << " SOS programming tests passed.\n";

    cout << "--------------------------------------------------------------\n\n";
    return make_tuple(test_count - fail_count, test_count);
}

void run_tests() {
    cout << "Running in test mode...\n";
    cout << "--------------------------------------------------------------\n\n";

    tuple<int, int> parse_results, program_results_scs, program_results_mosek, program_results_sosadmm;
    parse_results = run_polynomial_parser_tests();

//    string solver_choice = "scs";
//    program_results_scs = run_sos_program_tests(solver_choice);
//    solver_choice = "mosek";
//    program_results_mosek = run_sos_program_tests(solver_choice);
    string solver_choice = "sosadmm";
    program_results_sosadmm = run_sos_program_tests(solver_choice);

    cout << "Test summary:\t" << get<0>(parse_results) << "/" << get<1>(parse_results) << " parsing tests passed.\n";
//    cout << "\t\t" << get<0>(program_results_scs) << "/" << get<1>(program_results_scs) << " SOS programming tests passed with SCS.\n";
//    cout << "\t\t" << get<0>(program_results_mosek) << "/" << get<1>(program_results_mosek) << " SOS programming tests passed with MOSEK.\n";
    cout << "\t\t" << get<0>(program_results_sosadmm) << "/" << get<1>(program_results_sosadmm) << " SOS programming tests passed with MOSEK.\n";
    cout << "--------------------------------------------------------------\n";
}
