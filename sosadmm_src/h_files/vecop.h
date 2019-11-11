//
// Created by Jan Schappi on 02.11.19.
//

#ifndef SOSADMM_BOLT_VECOP_H
#define SOSADMM_BOLT_VECOP_H

#include <vector>
#include <thread>
#include <iostream>
#include <cmath>
using namespace std;

double norm(vector<double>& v);
double normsq(vector<double>& v);
double maximum(double d1, double d2, double d3, double d4);
double normOfDiffsq(vector<double>& v1,vector<double>& v2);

void addAndMulSecond(vector<double>& target, vector<double>& v1, vector<double>& v2, double mul, int NTHREADS);
void addAndMulSecond_thr(vector<double>& target, vector<double>& v1, vector<double>& v2, double mul, int first, int lastPlus);
void addAndMulSecond_par(vector<double>& target, vector<double>& v1, vector<double>& v2, double mul, int NTHREADS);
void addAndMulSecond_seq(vector<double>& target, vector<double>& v1, vector<double>& v2, double mul);

void sub(vector<double>& target, vector<double>& v1, vector<double>& v2, int NTHREADS);
void sub_thr(vector<double>& target, vector<double>& v1, vector<double>& v2, int first, int lastPlus);
void sub_par(vector<double>& target, vector<double>& v1, vector<double>& v2, int NTHREADS);
void sub_seq(vector<double>& target, vector<double>& v1, vector<double>& v2);

void mulPointwise(vector<double>& target, vector<double>& v1, vector<double>& v2, int NTHREADS);
void mulPointwise_thr(vector<double>& target, vector<double>& v1, vector<double>& v2, int first, int lastPlus);
void mulPointwise_par(vector<double>& target, vector<double>& v1, vector<double>& v2, int NTHREADS);
void mulPointwise_seq(vector<double>& target, vector<double>& v1, vector<double>& v2);

void addAndDivSecond_trivial(vector<double>& target, vector<double>& v1, vector<double>& v2, double div, int first, int lastPlus);
void addAndDivSecond_iterator(vector<double>& target, vector<double>& v1, vector<double>& v2, double div, int first, int lastPlus);

template <typename T>
double avg(vector<T>& v){
    T* pos = &v[0];
    T* lastPlus = &v[v.size()];
    double sum = 0;
    while(pos < lastPlus){
        sum += *(pos);
        ++pos;
    }
    return sum/v.size();
}

template <typename  T>
void abs(vector<T>& v){
    T* pos = &v[0];
    T* lastPlus = &v[v.size()];
    double sum = 0;
    while(pos < lastPlus){
        if (*(pos) < 0) *(pos) = - *(pos);
        ++pos;
    }
}

template <typename T>
void printVectorNoNewLine(vector<T>& v, string event){
    cout << event << " [";
    for(int i = 0; i < v.size(); ++i){
        cout << v[i] << " ";
    }
    cout << "]";
}

#endif //SOSADMM_BOLT_VECOP_H
