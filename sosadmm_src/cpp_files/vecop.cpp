//
// Created by Jan Schappi on 02.11.19.
//
#include "vecop.h"



// ADDANDMULSECOND
void addAndMulSecond(vector<double>& target, vector<double>& v1, vector<double>& v2, double mul, int NTHREADS){
    if(target.size() != v1.size() || target.size() != v2.size()){
        cout << "Error in addAndDivSecond: Vectors must be of equal size." << target.size() << v1.size() << v2.size() << endl;
        exit(1);
    }
    if(NTHREADS <= 1) addAndMulSecond_seq(target,v1,v2,mul);
    else addAndMulSecond_par(target,v1,v2,mul,NTHREADS);
}


void addAndMulSecond_par(vector<double>& target, vector<double>& v1, vector<double>& v2, double mul, int NTHREADS){
    int blocksize = target.size()/NTHREADS;
    int first,last;
    vector<thread> threads;
    for(int th = 0; th < NTHREADS-1; ++th){
        first = th*blocksize; last = (th+1)*blocksize;
        threads.push_back(thread(addAndMulSecond_thr,ref(target),ref(v1),ref(v2),mul,first,last));
    }
    threads.push_back(thread(addAndMulSecond_thr,ref(target),ref(v1),ref(v2),mul,last,target.size()));
    for(int th = 0; th < NTHREADS; ++th){
        threads[th].join();
    }
}

void addAndMulSecond_thr(vector<double>& target, vector<double>& v1, vector<double>& v2, double mul, int first, int lastPlus){
    double* pos1 = &v1[first]-1;
    double* pos2 = &v2[first]-1;
    double* postar = &target[first]-1;
    double* lastPlus_tar = &target[lastPlus]-1;
    while(postar < lastPlus_tar){
        *(++postar) = *(++pos1) + *(++pos2)*mul;
    }
}

void addAndMulSecond_seq(vector<double>& target, vector<double>& v1, vector<double>& v2, double mul){
    double* pos1 = &v1[0]-1;
    double* pos2 = &v2[0]-1;
    double* postar = &target[0]-1;
    double* lastPlus_tar = &target[target.size()]-1;
    while(postar < lastPlus_tar){
        *(++postar) = *(++pos1) + *(++pos2)*mul;
    }
}

// SUBTRACTION
void sub(vector<double>& target, vector<double>& v1, vector<double>& v2, int NTHREADS){
    if(target.size() != v1.size() || target.size() != v2.size()){
        cout << "Error in sub: Vectors must be of equal size." << target.size() << v1.size() << v2.size() << endl;
        exit(1);
    }
    if(NTHREADS <= 1) sub_seq(target,v1,v2);
    else sub_par(target,v1,v2,NTHREADS);
}


void sub_seq(vector<double>& target, vector<double>& v1, vector<double>& v2){
    double* pos1 = &v1[0]-1;
    double* pos2 = &v2[0]-1;
    double* postar = &target[0]-1;
    double* lastPlus_tar = &target[target.size()]-1;
    while(postar < lastPlus_tar){
        *(++postar) = *(++pos1) - *(++pos2);
    }
}

void sub_par(vector<double>& target, vector<double>& v1, vector<double>& v2, int NTHREADS){
    int blocksize = target.size()/NTHREADS;
    int first,last;
    vector<thread> threads;
    for(int th = 0; th < NTHREADS-1; ++th){
        first = th*blocksize; last = (th+1)*blocksize;
        threads.push_back(thread(sub_thr,ref(target),ref(v1),ref(v2),first,last));
    }
    threads.push_back(thread(sub_thr,ref(target),ref(v1),ref(v2),last,target.size()));
    for(int th = 0; th < NTHREADS; ++th){
        threads[th].join();
    }
}

void sub_thr(vector<double>& target, vector<double>& v1, vector<double>& v2, int first, int lastPlus){
    double* pos1 = &v1[first]-1;
    double* pos2 = &v2[first]-1;
    double* postar = &target[first]-1;
    double* lastPlus_tar = &target[lastPlus]-1;
    while(postar < lastPlus_tar){
        *(++postar) = *(++pos1) - *(++pos2);
    }
}

//POINTWISE MULTIPLICATION
void mulPointwise(vector<double>& target, vector<double>& v1, vector<double>& v2, int NTHREADS){
    if(target.size() != v1.size() || target.size() != v2.size()){
        cout << "Error in mulPointwise: Vectors must be of equal size." << target.size() << v1.size() << v2.size() << endl;
        exit(1);
    }
    if(NTHREADS <= 1) mulPointwise_seq(target,v1,v2);
    else mulPointwise_par(target,v1,v2,NTHREADS);
}


void mulPointwise_par(vector<double>& target, vector<double>& v1, vector<double>& v2, int NTHREADS){
    int blocksize = target.size()/NTHREADS;
    int first,last;
    vector<thread> threads;
    for(int th = 0; th < NTHREADS-1; ++th){
        first = th*blocksize; last = (th+1)*blocksize;
        threads.push_back(thread(mulPointwise_thr,ref(target),ref(v1),ref(v2),first,last));
    }
    threads.push_back(thread(mulPointwise_thr,ref(target),ref(v1),ref(v2),last,target.size()));
    for(int th = 0; th < NTHREADS; ++th){
        threads[th].join();
    }
}

void mulPointwise_thr(vector<double>& target, vector<double>& v1, vector<double>& v2, int first, int lastPlus){
    double* pos1 = &v1[first]-1;
    double* pos2 = &v2[first]-1;
    double* postar = &target[first]-1;
    double* lastPlus_tar = &target[lastPlus]-1;
    while(postar < lastPlus_tar){
        *(++postar) = *(++pos1) * *(++pos2);
    }
}

void mulPointwise_seq(vector<double>& target, vector<double>& v1, vector<double>& v2){
    double* pos1 = &v1[0]-1;
    double* pos2 = &v2[0]-1;
    double* postar = &target[0]-1;
    double* lastPlus_tar = &target[target.size()]-1;
    while(postar < lastPlus_tar){
        *(++postar) = *(++pos1) * *(++pos2);
    }
}

// USELESS
void addAndDivSecond_trivial(vector<double>& target, vector<double>& v1, vector<double>& v2, double div, int first, int lastPlus){
    for(int i = 0; i < target.size(); ++i){
        target[i] = v1[i]+v2[i]/div;
    }
}

void addAndDivSecond_iterator(vector<double>& target, vector<double>& v1, vector<double>& v2, double div, int first, int lastPlus){
    std::transform(v1.begin()+first,v1.begin()+lastPlus,v2.begin(),target.begin()+first,plus<double>());
}

double normsq(vector<double>& v){
    double* pos = &v[0];
    double* lastPlus = &v[v.size()];
    double sum = 0;
    while(pos < lastPlus){
        sum += *(pos)* *(pos);
        ++pos;
    }
    return sum;
}

double norm(vector<double>& v){
    return sqrt(normsq(v));
}

double maximum(double d1, double d2, double d3, double d4){
    return max(max(d1,d2),max(d3,d4));
}

double normOfDiffsq(vector<double>& v1,vector<double>& v2){
    if(v1.size() != v2.size()){
        cout << "Error in norm of Diff: vectors must be of same size!" << endl;
        exit(1);
    }
    vector<double> v3(v1.size());
    sub(v3,v1,v2,1);
    return normsq(v3);
}