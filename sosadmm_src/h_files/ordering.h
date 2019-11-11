//
// Created by Jan Schappi on 02.11.19.
//

#ifndef SOSADMM_BOLT_ORDERING_H
#define SOSADMM_BOLT_ORDERING_H

#include <vector>
#include <iostream>
using namespace std;

int positionOfAlpha(const vector<int>& alpha);
int sum(const vector<int>& alpha,const int start,const int end);
unsigned long int n_choose_k(int a, int b);
unsigned long int nMonomials(int N, int d);
vector<vector<int>> generate_all_exponents(int N, int d);

template<typename T>
void addSmall(vector<T>& result, const vector<T>& a, const vector<T>& b);

template<typename T>
void printVector(vector<T> v, string name){
    cout << name << ": [";
    for(int i = 0; i < v.size(); ++i){
        cout << v[i] << " ";
    }
    cout << "]" << endl;
}

template<typename T>
void addSmall(vector<T>& result, const vector<T>& a, const vector<T>& b){
    if(a.size() != b.size() || a.size() != result.size()){
        cout << "Error in sum vectors: Dimensions must be the same." << endl;
        exit(1);
    }
    for(int i = 0; i<result.size(); ++i){
        result[i] = a[i]+b[i];
    }
}

#endif //SOSADMM_BOLT_ORDERING_H
