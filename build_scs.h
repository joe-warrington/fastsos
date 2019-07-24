//
// Created by Joe Warrington on 2019-07-23.
//

#ifndef FASTSOS_BUILD_SCS_H
#define FASTSOS_BUILD_SCS_H

#endif //FASTSOS_BUILD_SCS_H

static const char *simple_lp(void);

int size_of_vmat(int side_length);

tuple<ScsCone, ScsData, ScsSolution, ScsInfo> create_scs_model(PolyInfo& f_info, vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos,
                                                      int dim, int deg, int total_vars, string positivity_condition);

void create_scs_coeff_matches(tuple<ScsCone, ScsData, ScsSolution, ScsInfo> &M,
                              vector<double> &f_mono_coeffs, vector<vector<int> > &f_mono_exponents,
                              vector<vector<double> > &g_mono_coeffs, vector<vector<vector<int> > > &g_mono_exponents,
                              vector<vector<double> > &h_mono_coeffs, vector<vector<vector<int> > > &h_mono_exponents,
                              int dim, int deg, unsigned long int s_of_d, vector<unsigned long int> &s_of_d_minus_djs,
                              vector<unsigned long int> &s_of_d_minus_dj2s, int output_level);

tuple<int, vector<int>, vector<string> > calc_n_vars(vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos, int n, int d);

int compute_scs_nonzeros(PolyInfo& f_info, vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos,
        int dim, int deg, vector<vector<int> > exp_list_2d);