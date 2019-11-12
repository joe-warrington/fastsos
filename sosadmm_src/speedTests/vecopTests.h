//
// Created by Jan Schappi on 02.11.19.
//

#ifndef SOSADMM_BOLT_VECOPTESTS_H
#define SOSADMM_BOLT_VECOPTESTS_H

#include "vecop.h"
#include "timing.h"
#include "ordering.h"

void checkIfEquals(vector<double> test, double equal){
    for(int i = 0; i < test.size(); ++i){
        if(test[i] != equal){
            cout << "Vector is not equal to " << equal << " at position " << i << endl;
            return;
        }
    }
    cout << "Vector is equal to " << equal << endl;
}

void testSumAndDiv(){
    const int length = 10000000;
    vector<double> v1(length);
    vector<double> v2(length);
    vector<double> target(length);
    for(int i = 0; i < length; ++i){
        v1[i] = 2;
        v2[i] = 3;
    }
    timestamp_t t0 = timenow();
    addAndMulSecond(target,v1,v2,2,1);
    timestamp_t t1 = timenow();
    cout << "Single Thread SumAndMult " << time_string(t1-t0) << endl;
    checkIfEquals(target,8);
    t0 = timenow();
    addAndMulSecond(target,v1,v2,2,2);
    t1 = timenow();
    cout << "Two Thread SumAndMult " << time_string(t1-t0) << endl;
    checkIfEquals(target,8);

    t0 = timenow();
    sub(target,v1,v2,1);
    t1 = timenow();
    cout << "Single Thread sub " << time_string(t1-t0) << endl;
    checkIfEquals(target,-1);

    t0 = timenow();
    sub(target,v1,v2,2);
    t1 = timenow();
    cout << "Two Thread sub " << time_string(t1-t0);
    checkIfEquals(target,-1);

    t0 = timenow();
    mulPointwise(target,v1,v2,1);
    t1 = timenow();
    cout << "Single Thread sub " << time_string(t1-t0) << endl;
    checkIfEquals(target,6);

    t0 = timenow();
    mulPointwise(target,v1,v2,2);
    t1 = timenow();
    cout << "Two Thread sub " << time_string(t1-t0);
    checkIfEquals(target,6);

}

#endif //SOSADMM_BOLT_VECOPTESTS_H
