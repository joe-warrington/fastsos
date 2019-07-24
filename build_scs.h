//
// Created by Joe Warrington on 2019-07-23.
//

#ifndef FASTSOS_BUILD_SCS_H
#define FASTSOS_BUILD_SCS_H

#endif //FASTSOS_BUILD_SCS_H

static const char *simple_lp(void);

int size_of_vmat(int side_length);

tuple<int, vector<int>, vector<string> > calc_n_vars(vector<PolyInfo>& g_infos, vector<PolyInfo>& h_infos, int n, int d);