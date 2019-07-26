// (c) 2019 ETH Zurich, Automatic Control Lab, Joe Warrington

#include "mosek.h"
#include "fusion.h"
using namespace mosek::fusion;

#ifndef FASTSOS_SOS_H
#define FASTSOS_SOS_H

#endif //FASTSOS_SOS_H

unsigned long int a_choose_b(int a, int b);

unsigned long int n_monomials(int n, int d);

unsigned long int get_deg_start(int n, int deg_this_level);

template <class T>
void print_vec(const string& string_in, const T& v);

string monomial_to_string(const vector<int> v_in);

vector<int> add_vecs(const vector <int> v1, const vector <int> v2);

vector<int> sub_vecs(const vector <int> v1, const vector <int> v2);

template <class T>
T sum_vec(const vector <T>& v);

void eliminate_unused_dims(int& n, vector<vector<int> >& f_exps, vector<vector<vector<int> > >& g_exps_list,
                           vector<vector<vector<int> > >& h_exps_list, int output_level);

vector<vector <int> > generate_all_exponents(int n, int d, int output_level);

int compute_legal_d(PolyInfo f_info, vector<PolyInfo> g_infos, vector<PolyInfo> h_infos, int d_request);

tuple<double, int, string> sos_level_d(string& f_string, vector<string>& g_strings, vector<string>& h_strings,
                               int d_request, string& positivity_condition, int output_level, string solver_choice);