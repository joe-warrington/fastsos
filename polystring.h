//
// Created by Joe Warrington on 2019-05-31.
//

#ifndef GIT_CPPSOS_POLYSTRING_H
#define GIT_CPPSOS_POLYSTRING_H

#endif //GIT_CPPSOS_POLYSTRING_H

using namespace std;

struct PolyInfo{
    int dimension = 0;
    int degree = 0;
    int n_terms = 0;
};

PolyInfo infer_dim_and_n_terms(const std::string& s);

void parse_poly(string& s, vector<double>& mono_coeffs_in, vector<vector<int> >& mono_exponents_in,
                PolyInfo& info);