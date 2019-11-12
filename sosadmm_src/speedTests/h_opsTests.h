//
// Created by Jan Schappi on 02.11.19.
//

#ifndef SOSADMM_BOLT_H_OPSTESTS_H
#define SOSADMM_BOLT_H_OPSTESTS_H
//#include "sdp.h"
#include "h_ops.h"
#include <iostream>
using namespace std;

void makeIncr(vector<double>& x){
    for(int i = 0; i < x.size(); ++i){
        x[i] = i+1;
    }
}

void testing_h_ops(sdp& SDP){
    cout << "Testing h_ops." << endl;
    cout << "Make sure the input problem is ../problems/test_h_ops.txt" << endl;
    cout << "Should be [1,2,4,7,5,10,12,8,13,14,15,14,].";

    vector<double> x_sized(SDP.n,3);
    vector<double> z_sized(SDP.NNZ,4);
    vector<double> w_sized(SDP.m,5);

    makeIncr(z_sized);
    sumHiT(x_sized, z_sized, SDP.colsA);
    printVector(x_sized,"sumHiT(1 2 3 ..)");

    cout << "Should be "
            "1 2 12 \n"
            "3 5 11 \n"
            "4 8 11 \n"
            "6 12 \n"
            "7 9 \n"
            "10" << endl;
    makeIncr(x_sized);
    Hix_seq(z_sized,x_sized,SDP.colsA);
    printVector(z_sized,"z_sized");

    makeIncr(w_sized);
    Hiai_seq(z_sized, w_sized, SDP.rowPointA, SDP.valsA);
    cout << "Should be 1 1 -1 \n"
            "2 2 2 \n"
            "3 3 9 \n"
            "4 4 \n"
            "5 5 \n"
            "6 " << endl;
    printVector(z_sized,"z_sized");
    makeIncr(z_sized);
    HiaiT_seq(w_sized,z_sized,SDP.rowPointA,SDP.valsA);
    cout << "Should be [0,15,42,21,25,14,]" << endl;
    printVector(w_sized,"w");
}

#endif //SOSADMM_BOLT_H_OPSTESTS_H
