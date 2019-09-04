// (c) 2019 ETH Zurich, Automatic Control Lab, Joe Warrington

#include "util.h"  // includes cones.h and scs.h; scs.h includes glbopts.h and aa.h
#include "amatrix.h"  // includes glbopts.h

#ifndef FASTSOS_BUILD_SCS_H
#define FASTSOS_BUILD_SCS_H

#endif //FASTSOS_BUILD_SCS_H

int size_of_vmat(int side_length);

tuple<ScsCone *, ScsData *, ScsSolution *, ScsInfo, int> create_scs_model(
        PolyInfo& f_info, vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos,
        vector<vector<vector<int> > > &g_mono_exponents, vector<vector<vector<int> > > &h_mono_exponents,
        int dim, int deg, vector<vector<int> > exp_list_2d,
        int total_vars, vector<unsigned long int> &s_of_d_minus_djs,
        vector<unsigned long int> &s_of_d_minus_dj2s, string positivity_condition, int output_level);

void create_scs_coeff_matches(tuple<ScsCone *, ScsData *, ScsSolution *, ScsInfo, int> M,
                              vector<double> &f_mono_coeffs, vector<vector<int> > &f_mono_exponents,
                              vector<vector<double> > &g_mono_coeffs, vector<vector<vector<int> > > &g_mono_exponents,
                              vector<vector<double> > &h_mono_coeffs, vector<vector<vector<int> > > &h_mono_exponents,
                              int dim, int deg, vector<vector<int> > exp_list_2d,
                              unsigned long int s_of_d, vector<unsigned long int> &s_of_d_minus_djs,
                              vector<unsigned long int> &s_of_d_minus_dj2s, int output_level);

tuple<int, vector<int>, vector<string> > calc_n_vars(vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos, int n, int d);

int compute_scs_nonzeros(PolyInfo& f_info, vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos,
                         vector<vector<vector<int> > > & g_mono_exponents, vector<vector<vector<int> > > & h_mono_exponents,
                         int dim, int deg, const vector<vector<int> > & exp_list_2d);

int compute_scs_nonzeros_new(PolyInfo& f_info, vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos,
                             vector<vector<vector<int> > > & g_mono_exponents, vector<vector<vector<int> > > & h_mono_exponents,
                             int dim, int deg, const vector<vector<int> > & exp_list_2d, vector<unsigned long int> &s_of_d_minus_djs,
                             vector<unsigned long int> &s_of_d_minus_dj2s);

tuple<double, int, string> solve_with_scs(tuple<int, int, int, PolyInfo, vector<PolyInfo>, vector<PolyInfo>,
        vector<double>, vector<vector<int> >,
        vector<vector <double> >, vector<vector <vector <int> > >,
        vector<vector <double> >, vector<vector <vector <int> > > > data_tuple, int d, string positivity_condition, int output_level);