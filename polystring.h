// (c) 2019 ETH Zurich, Automatic Control Lab, Joe Warrington

#ifndef FASTSOS_POLYSTRING_H
#define FASTSOS_POLYSTRING_H

#endif //FASTSOS_POLYSTRING_H

using namespace std;

struct PolyInfo{
    int dimension = 0;
    int degree = 0;
    int n_terms = 0;
};

PolyInfo infer_dim_and_n_terms(const std::string& s, bool diag_msg);

void print_polynomial_table(vector<double>& mono_coeffs, vector<vector<int> >& mono_exponents);

void parse_poly(string& s, vector<double>& mono_coeffs_in, vector<vector<int> >& mono_exponents_in,
                PolyInfo& info, bool diag_msg);

tuple<string, vector<string>, vector<string>, bool > read_problem_from_file(const string& filename);