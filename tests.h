//
// Created by Joe Warrington on 2019-06-27.
//

#ifndef GIT_CPPSOS_TESTS_H
#define GIT_CPPSOS_TESTS_H

#endif //GIT_CPPSOS_TESTS_H

tuple<bool, bool, bool, PolyInfo, vector<double>, vector<vector<int> >> test_polynomial_parse(
        string s, vector<int> info_in, vector<double> coeffs_in, vector<vector<int> > exponents_in);

tuple<int, int> run_polynomial_parser_tests();

tuple<bool, bool, bool, bool, double> test_sos_program(string& obj_in, vector<string> constrs_in, int d_in, string pos_cert_in,
                                                       tuple<bool, bool, double, double, bool> expected_result_in);

tuple<int, int> run_sos_program_tests();

void run_tests();