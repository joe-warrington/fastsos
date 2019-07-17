//
// Created by Joe Warrington on 2019-06-01.
//

#include "mosek.h"
#include "fusion.h"
using namespace mosek::fusion;

#ifndef GIT_CPPSOS_SOS_H
#define GIT_CPPSOS_SOS_H

#endif //GIT_CPPSOS_SOS_H

string time_string(unsigned long long us_in);

unsigned long int a_choose_b(int a, int b);

unsigned long int n_monomials(int n, int d);

void eliminate_unused_dims(int& n, vector<vector<int> >& f_exps, vector<vector<vector<int> > >& g_exps_list, int output_level);

vector<vector <int> > generate_all_exponents(int n, int d, int output_level);

void constrain_to_cone(Model::t& M, Variable::t& matrix_var, int matrix_size, string& cone_type);

tuple<Model::t, Variable::t, Variable::t, vector<Variable::t>, vector<unsigned long int> > create_mosek_model(
        PolyInfo& f_info, vector<PolyInfo>& g_infos, int n, int d, unsigned long int s_of_d, string positivity_condition);

void create_coeff_matches(Model::t& M, vector<double>& f_mono_coeffs, vector<vector<int> >& f_mono_exponents,
                          vector<vector <double> >& g_mono_coeffs, vector<vector <vector <int> > >& g_mono_exponents,
                          int n, int d, unsigned long int s_of_d, Variable::t& lambda, Variable::t& sigma_0,
                          vector<Variable::t>& sigma_j, vector<unsigned long int>& s_of_d_minus_djs, int output_level);

int compute_d(PolyInfo f_info, vector<PolyInfo> g_infos, int d_request);

template <class T>
void print_vec(const string& string_in, const T& v);

vector<int> add_vecs(vector <int> v1, vector <int> v2);
vector<int> sub_vecs(vector <int> v1, vector <int> v2);

template <class T>
T sum_vec(const vector <T>& v);

tuple<double, ProblemStatus, SolutionStatus, SolutionStatus> sos_level_d(
        string& f_string, vector<string>& g_strings, int d_request, string& positivity_condition, int output_level);