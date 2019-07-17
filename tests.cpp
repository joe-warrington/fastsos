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
        cout << "WARNING: " << fail_count << "/" << test_count << " polynomial parsing tests failed." << endl;
    else
        cout << "All " << test_count << " polynomial parsing tests passed.\n";

    cout << "--------------------------------------------------------------\n\n";
    return make_tuple(test_count - fail_count, test_count);
}

tuple<bool, bool, bool, bool, double> test_sos_program(
        string& obj_in, vector<string> constrs_in, int d_in, string pos_cert_in,
        tuple<bool, bool, double, double, bool> expected_result_in) {

    bool obj_val_correct = false;
    bool primal_status_correct = false;
    bool dual_status_correct = false;
    bool problem_status_correct = false;

    auto sol_info = sos_level_d(obj_in, constrs_in, d_in, pos_cert_in, 0);

    double obj_val = get<0>(sol_info);
    ProblemStatus problem_status = get<1>(sol_info);
    SolutionStatus solution_status_primal = get<2>(sol_info);
    SolutionStatus solution_status_dual = get<3>(sol_info);

    bool expected_primal_feasible = get<0>(expected_result_in);
    bool expected_primal_bounded = get<1>(expected_result_in);
    double expected_obj_val = get<2>(expected_result_in);
    double expected_tolerance = get<3>(expected_result_in);
    bool expected_num_problems = get<4>(expected_result_in);

    bool report_obj_value = false;

    // Check overall status returned by problem and set correctness flags accordingly
    switch (problem_status) {
        case ProblemStatus::PrimalAndDualFeasible :
            report_obj_value = true;
            if (solution_status_primal == SolutionStatus::Optimal && solution_status_dual == SolutionStatus::Optimal) {
                if (abs(obj_val - expected_obj_val) <= expected_tolerance) {
                    obj_val_correct = true;
                }
            }
            if (expected_primal_feasible) primal_status_correct = true;
            if (expected_primal_bounded) dual_status_correct = true;
            if (expected_primal_feasible && expected_primal_bounded) problem_status_correct = true;
            break;
        case ProblemStatus::PrimalFeasible :
            if (solution_status_primal == SolutionStatus::Optimal && solution_status_dual == SolutionStatus::Optimal) {
                if (abs(obj_val - expected_obj_val) <= expected_tolerance) {
                    obj_val_correct = true;
                }
            }
            break;
        case ProblemStatus::PrimalInfeasible :
            obj_val_correct = true; // Objective value deemed correct if problem is infeasible or unbounded
            break;
        case ProblemStatus::PrimalInfeasibleOrUnbounded :
            break;
        case ProblemStatus::DualFeasible :
            break;
        case ProblemStatus::DualInfeasible :
            if (expected_primal_feasible) primal_status_correct = true;
            if (!expected_primal_bounded) dual_status_correct = true;
            if (expected_primal_feasible && !expected_primal_bounded) problem_status_correct = true;
            obj_val_correct = true; // Objective value deemed correct if problem is infeasible or unbounded
            break;
        case ProblemStatus::PrimalAndDualInfeasible :
            obj_val_correct = true; // Objective value deemed correct if problem is infeasible or unbounded
            break;
        case ProblemStatus::IllPosed :
            cout << "";
            break;
        case ProblemStatus::Unknown :
            if (expected_num_problems) problem_status_correct = true;
            primal_status_correct = true;
            dual_status_correct = true;
            obj_val_correct = true; // Objective value deemed correct if problem is infeasible or unbounded
            break;
    }

    if (report_obj_value) cout << "Objective value:\t" << obj_val << endl;
    cout << "Problem status:\t\t" << problem_status << endl;
    cout << "Primal solution status:\t" << solution_status_primal << endl;
    cout << "Dual solution status:\t" << solution_status_dual << endl;

    return make_tuple(obj_val_correct, primal_status_correct, dual_status_correct, problem_status_correct, obj_val);
}

tuple<int, int> run_sos_program_tests() {
    cout << "Running SOS program tests..." << endl;
    vector<string> input_objs = {"x1",
                                 "x1",
                                 "x1",
                                 "0.5 + x1^1",
                                 "x2",
                                 "x1 + x2"};
    vector<vector<string> > input_constrs = {{"x1 + 0.5", "0.5 - x1"},
                                             {"-x1"},
                                             {"-x1", "x1 - 1"},
                                             {"-x1^2 - x2^2 + 1"},
                                             {"-x2^2 - x5^2 + 1"},
                                             {"-x1^2 - x2^2 + 1"}};
    vector<int> relaxation_degree = {2,
                                     2,
                                     2,
                                     2,
                                     2,
                                     2};
    vector<string> positivity_certs = {"PSD",
                                       "PSD",
                                       "PSD",
                                       "PSD",
                                       "PSD",
                                       "PSD"};
    vector<bool> primal_feasible = {true,
                                    true,
                                    false,
                                    true,
                                    true,
                                    true};
    vector<bool> primal_bounded = {true,
                                   false,
                                   true,
                                   true,
                                   true,
                                   true};
    vector<double> opt_values = {-0.5,
                                 0.0,
                                 0.0,
                                 -0.5,
                                 -1.0,
                                 -1.41421};
    vector<double> sol_tols = {1e-6,
                               1e-6,
                               1e-6,
                               1e-4,
                               1e-4,
                               1e-4};
    vector<bool> num_problems = {false,
                                 true,
                                 false,
                                 false,
                                 false,
                                 false};

    int n_tests = input_objs.size();

    int test_count = 0;
    int fail_count = 0;
    tuple<bool, bool, double, double, bool> expected_result;
    tuple<bool, bool, bool, bool, double> test_output;

    for (int i = 0; i < n_tests; i++) {
        test_count++;
        string problem_string = "min " + input_objs[i];
        if (input_constrs[i].size() > 0) problem_string += "   s.t.  ";
        for (int j = 0; j < input_constrs[i].size(); j++)
            problem_string += input_constrs[i][j] + " >= 0,  ";
        cout << "#" << i + 1 << "\t" << problem_string << "\b\b" << endl;
        expected_result = make_tuple(primal_feasible[i], primal_bounded[i], opt_values[i], sol_tols[i], num_problems[i]);

        test_output = test_sos_program(input_objs[i], input_constrs[i], relaxation_degree[i], positivity_certs[i], expected_result);

        if (!get<0>(test_output) || !get<1>(test_output) || !get<2>(test_output) || !get<3>(test_output)) {
            cout << color_red << "FAIL: ";
            if (!get<0>(test_output)) {
                cout << "Objective value is incorrect: " << get<4>(test_output) << " instead of "
                        << opt_values[i] << " +- " << sol_tols[i] << endl;
            }
            if (!get<1>(test_output)) {
                cout << "Primal feasibility status is incorrect." << endl;
            }
            if (!get<2>(test_output)) {
                cout << "Primal boundedness status is incorrect." << endl;
            }
            if (!get<3>(test_output)) {
                cout << "Overall solution status is incorrect." << endl;
            }
            fail_count++;
            cout << color_reset << "\n";
        } else
            cout << color_green << "PASS: solution behaved as expected.\n" << color_reset << endl;
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
    tuple<int, int> parse_results, program_results;
    parse_results = run_polynomial_parser_tests();
    program_results = run_sos_program_tests();
    cout << "Test summary:\t" << get<0>(parse_results) << "/" << get<1>(parse_results) << " parsing tests passed.\n";
    cout << "\t\t" << get<0>(program_results) << "/" << get<1>(program_results) << " SOS programming tests passed.\n";
    cout << "--------------------------------------------------------------\n";
}
