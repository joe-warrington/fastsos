//
// Created by Jan Schappi on 03.11.19.
//

#ifndef SOSADMM_BOLT_PROJECTIONTESTS_H
#define SOSADMM_BOLT_PROJECTIONTESTS_H
#include <vector>
#include "timing.h"
using namespace std;

void addvvT2(vector<double>& x, double* z, double* w, int m, int N, int i){
    double* posI = &z[N*i]-1; double* posILastPlus = posI+N;
    double* posJ = posI;
    double* posX = &x[0]-1;
    double eig = w[i];
    int skip = 0;
    while(posI < posILastPlus){
        posJ = posI;
        ++posI;
        while(posJ < posILastPlus){
            *(++posX) += *(posI) * *(++posJ) * eig;
        }
        posX += (++skip);
    }
}

void makeSymFromHalf2(vector<double>& x, int N){
    for(int i = 0; i < N; ++i){
        for(int j = i+1; j < N; ++j){
            x[i+N*j] = x[j + N*i];
        }
    }

}

void writeBack2(vector<double>& x, double* z, double* w, int m, int N){
    memset(&x[0],0,x.size()*sizeof(double));
    for(int i = 0; i < m; ++i){
        addvvT2(x, z, w, m, N, i);
    }
    makeSymFromHalf2(x,N);
}

void testWriteBack(){
    int N = 231;
    double w[N];
    for(int i = 0; i < N; ++i){
        w[i] = 3;
    }
    double z[N*N];
    for(int i = 1; i < N; ++i){
        for(int j = 1; j < N; ++j){
            w[i] = i-j;
        }
    }
    vector<double> x(N*N,0);
    timestamp_t t0 = timenow();
    writeBack2(x,z,w,N,N);
    cout << time_string(timenow()-t0);
}

void testPSD(){

    vector<double> y(9);
    for(int i = 0; i < 9; ++i) y[i] = i+1;
    printVector(y,"iVector");
    projectOnCone_LAPACK(y,0,3);
    printVector(y,"After");

    int N = 231;
    vector<double> x(N*N);

    std::random_device rd{};
    std::mt19937 gen{rd()};

    //std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,1);

    for(int test = 0; test < 10; ++test){
        for(int i = 0; i < N; ++i){
            for(int j = i; j < N; ++j){
                x[N*i+j] = distribution(gen);
                x[N*j+i] = x[N*j+i];
            }
        }
        timestamp_t t0 = timenow();
        projectOnCone_LAPACK(x,0,N);
        t0 = timenow()-t0;
        cout << time_string(t0) << endl;
    }



}

#endif //SOSADMM_BOLT_PROJECTIONTESTS_H
