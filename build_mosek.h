//
// Created by Joe Warrington on 2019-07-23.
//

#include "mosek.h"
#include "fusion.h"
using namespace mosek::fusion;

#ifndef FASTSOS_BUILD_MOSEK_H
#define FASTSOS_BUILD_MOSEK_H

#endif //FASTSOS_BUILD_MOSEK_H


void mosek_constrain_to_cone(Model::t &M, Variable::t &matrix_var, int matrix_size, const string &cone_type);

tuple<Model::t, Variable::t, Variable::t, vector<Variable::t>, vector<unsigned long int>,
        vector<Variable::t>, vector<unsigned long int> > create_mosek_model(
        PolyInfo& f_info, vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos,
        int n, int d, unsigned long int s_of_d, string positivity_condition);

void create_mosek_coeff_matches(Model::t &M, vector<double> &f_mono_coeffs, vector<vector<int> > &f_mono_exponents,
                                vector<vector<double> > &g_mono_coeffs, vector<vector<vector<int> > > &g_mono_exponents,
                                vector<vector<double> > &h_mono_coeffs, vector<vector<vector<int> > > &h_mono_exponents,
                                int n, int d, unsigned long int s_of_d, Variable::t &lambda, Variable::t &sigma_0,
                                vector<Variable::t> &sigma_j, vector<unsigned long int> &s_of_d_minus_djs,
                                vector<Variable::t> &tau_j, vector<unsigned long int> &s_of_d_minus_dj2s,
                                int output_level);

tuple<double, int, string> solve_with_mosek(tuple<int, int, int, PolyInfo, vector<PolyInfo>, vector<PolyInfo>,
                                               vector<double>, vector<vector<int> >,
                                               vector<vector <double> >, vector<vector <vector <int> > >,
                                               vector<vector <double> >, vector<vector <vector <int> > > > data_tuple,
                                               int d, string positivity_condition, int output_level);
