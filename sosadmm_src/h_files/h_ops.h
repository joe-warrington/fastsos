//
// Created by Jan Schappi on 02.11.19.
//

#ifndef SOSADMM_BOLT_H_OPS_H
#define SOSADMM_BOLT_H_OPS_H

#include "sdp.h"
#include <thread>

void sumHiT(vector<double>& target_x, const vector<double>& arg_z, const vector<int>& colsA);
void Hix_seq(vector<double>& target_z, const vector<double>& input_x, const vector<int>& colsA);
void Hiai_seq(vector<double>& target_z, const vector<double>& input_w, const vector<int>& rowPointA, const vector<double>& valsA);
void HiaiT_seq(vector<double>& target_w, const vector<double>& input_z, const vector<int>& rowPointA, const vector<double>& valsA);
#endif //SOSADMM_BOLT_H_OPS_H
