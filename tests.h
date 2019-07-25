// (c) 2019 ETH Zurich, Automatic Control Lab, Joe Warrington

#ifndef FASTSOS_TESTS_H
#define FASTSOS_TESTS_H

#endif //FASTSOS_TESTS_H

tuple<bool, bool, bool, PolyInfo, vector<double>, vector<vector<int> >> test_polynomial_parse(
        string s, vector<int> info_in, vector<double> coeffs_in, vector<vector<int> > exponents_in);

tuple<int, int> run_polynomial_parser_tests();

tuple<bool, bool, double, int> test_sos_program(
        string& obj_in, vector<string> ineq_constrs_in, vector<string> eq_constrs_in, int d_in, string pos_cert_in,
        tuple<int, double, double, bool> expected_result_in, string solver_choice);

tuple<int, int> run_sos_program_tests(string solver_choice);

void run_tests();