//
// Created by Jan Schappi on 02.11.19.
//

#include "h_ops.h"

// This can NOT be done in parallel trivially due to summation
void sumHiT(vector<double>& target_x, const vector<double>& input_z, const vector<int>& colsA){

    memset(&target_x[0], 0, target_x.size() * sizeof target_x[0]);
    const int* posA = &colsA[0]-1;
    const double* posZ = &input_z[0]-1;
    const int* lastPlus_A = &colsA[colsA.size()]-1;
    while(posA < lastPlus_A){
        //--pos;
        //        //target_x[colsA[pos]] += input_z[pos];
        target_x[*(++posA)] += *(++posZ);
    }
}

// Just allocation - can be parallelized over colsA
void Hix_seq(vector<double>& target_z, const vector<double>& input_x, const vector<int>& colsA){
    const int* posA = &colsA[0]-1;
    const int* lastPlus_A = &colsA[colsA.size()]-1;
    double* posZ = &target_z[0]-1;
    while(posA < lastPlus_A){
        *(++posZ) = input_x[*(++posA)];
    }
}

// Just allocation - can be parallelized. Can be parallelized over i!
void Hiai_seq(vector<double>& target_z, const vector<double>& input_w, const vector<int>& rowPointA, const vector<double>& valsA){
    const double* posA = &valsA[0]-1;
    const double* lastPlus_A = &valsA[valsA.size()]-1;
    double* posZ = &target_z[0]-1;
    double wi = input_w[0];
    int i = 0;
    int counter = 0;
    const int* nextChange = &rowPointA[1]; // Always at least size 2 since it holds m+1 values
    while(posA < lastPlus_A){
        *(++posZ) = *(++posA) * wi;
        if(++counter == *(nextChange)){
            wi = input_w[++i];
            ++nextChange;
        }
    }
}

// Addition - not trivially parallelizable but probably possible (rowwise)
void HiaiT_seq(vector<double>& target_w, const vector<double>& input_z, const vector<int>& rowPointA, const vector<double>& valsA){
    memset(&target_w[0], 0, target_w.size() * sizeof target_w[0]);
    const double* posA = &valsA[0]-1;
    const double* lastPlus_A = &valsA[valsA.size()]-1;
    const double* posZ = &input_z[0]-1;
    double* posW = &target_w[0];
    int counter = 0;
    const int* nextChange = &rowPointA[1];
    while(posA < lastPlus_A){
        *(posW) += *(++posA)* *(++posZ);
        if(++counter == *(nextChange)){
            ++posW;
            ++nextChange;
        }
    }
}
