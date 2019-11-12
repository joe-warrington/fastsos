//
// Created by Jan Schappi on 03.11.19.
//

#ifndef SOSADMM_BOLT_UNCONSTRAINEDTEST_H
#define SOSADMM_BOLT_UNCONSTRAINEDTEST_H

#include <vector>
#include <random>
#include "vecop.h"
#include "ordering.h"
#include "admm.h"
using namespace std;

vector<double> squarePoly(vector<double>& coeff_d, vector<vector<int>>& alphas_d, int N, int d){
    int m = n_monomials(N,2*d);
    vector<double> coeff_2d(m,0);
    vector<int> alpha_lk(N);
    for(int l = 0; l < coeff_d.size(); ++l){
        for(int k = 0; k < coeff_d.size(); ++k){
            addSmall(alpha_lk,alphas_d[k],alphas_d[l]);
            coeff_2d[positionOfAlpha(alpha_lk)] += coeff_d[l]*coeff_d[k];
        }
    }
    return coeff_2d;
}

double norm(vector<double>& v, int start, int endPlus){
    double sum = 0;
    for(int i = start; i < endPlus; ++i){
        sum += v[i]*v[i];
    }
    return sqrt(sum);
}

// Returns coefficients for a polynomial with random normally distributed
vector<double> randn(int N, int d, double mean, double stddev){
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> sample(mean,stddev);
    int coeff_size = nMonomials(N,d);
    vector<double> coeffs(coeff_size);

    for(int i = 0; i < coeff_size; ++i){
        coeffs[i] = sample(gen);
    }
    return coeffs;
}

void randomPolynomialSolution(vector<double>& b, int N, int d){

    std::random_device rd{};
    std::mt19937 gen{rd()};

    //std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,1);


    vector<vector<int>> alpha = generate_all_exponents(N,2*d);

    int firstDeg_d = nMonomials(N,d-1);

    int monomials_N_d = nMonomials(N,d);
    vector<double> coeffsRest(monomials_N_d,0);
    vector<int> alpha_small(N);
    for(int n = 0; n < N; ++n){
        for(int l = 0; l < N; ++l) alpha_small[l] = 0;
        alpha_small[n] = d;
        coeffsRest[positionOfAlpha(alpha_small)] = 1;
    }
    b = squarePoly(coeffsRest,alpha,N,d);

    int entriesToFill = nMonomials(N,2*d-1);

    for(int i = 0; i < entriesToFill; ++i) b[i] = distribution(gen);
    double norm_b = norm(b,0,entriesToFill);
    for(int i = 0; i < entriesToFill; ++i) b[i] /= norm_b;
}

tuple<double,int,string,timestamp_t,double,double> callAdmm(PolyTable polytable, int d,double eps){
    timestamp_t t0 = timenow();
    sdp SDP(polytable,d);
    admm ADMM(SDP);
    ADMM.monitor = false;
    ADMM.eps = eps;
    tuple<double,int,string> sol = ADMM.min();
    t0 = timenow()-t0;
    //cout << get<2>(sol) << endl;
    return make_tuple(get<0>(sol),ADMM.itUsed,get<2>(sol),t0,ADMM.pres,ADMM.dres);
}

void run_unconstrained_experiment_with(int N, int d, int nExp, double eps){
    int dim_b = nMonomials(N,2*d);
    vector<double> b(dim_b,0);
    randomPolynomialSolution(b,N,d);
    vector<vector<int>> alpha = generate_all_exponents(N,2*d);
    PolyTable polytable;
    get<0>(polytable) = N;
    get<1>(polytable) = 0;
    get<2>(polytable) = 0;
    get<7>(polytable) = alpha;
    get<6>(polytable) = b;
    sdp SDP(polytable,d);
    admm ADMM(SDP);

    tuple<double,int,string,timestamp_t,double,double> eval;
    double avgErr = 0;

    for(int exp = 0; exp < nExp; ++exp){
        randomPolynomialSolution(b,N,d);
        get<6>(polytable) = b;
        eval = callAdmm(polytable,d,eps);


    }
}

void run_unconstrained_speed(){
    int d = 2;
    for(int N = 2; N <= 20; N += 2){
        run_unconstrained_experiment_with(N,d,1,1e-4);
    }
}

void run_PSD_projection(){
    std::random_device rd{};
    std::mt19937 gen{rd()};
    //std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,1);
    int N = 231;
    int l,k;
    vector<double> x(N*N);
    for(int i = 0; i < x.size(); ++i) {
        // i = kN+l j = lN+k
        x[i] = distribution(gen);
        l = i % N;
        k = i / N;
        x[k*N+l] = x[i];
    }
    timestamp_t t0;
    for(int i = 0; i < 100; ++i){
        t0 = timenow();
        projectOnCone_LAPACK(x,0,N);
        t0 = timenow()-t0;
    }

    cout << time_string(t0) << endl;
}

void totallyRandomPolynomial(vector<double>& coeffs, vector<vector<int>>& alphas, int N, int d){
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> sampler(0,1);
    for(int i = 0; i < coeffs.size(); ++i){
        coeffs[i] = sampler(gen);
    }
    alphas = generate_all_exponents(N,d);
}

void randomPolynomialSolution_SOS_Solution0(vector<double>& b, int N, int d){
    int nPolySOS = 4;
    int sz_d = nMonomials(N,d);
    int sz_2d = nMonomials(N,2*d);
    vector<double> polySq(sz_2d);
    vector<double> poly_d(sz_d);
    vector<vector<int>> alphas_d = generate_all_exponents(N,d);
    for(int i = 0; i < nPolySOS; ++i){
        poly_d = randn(N,d,0,1.0);
        poly_d[0] = 0;
        polySq = squarePoly(poly_d,alphas_d,N,d);
        addAndMulSecond(b,b,polySq,1,1);
    }
    //double norm_b = norm(b);
    //for(int i = 0; i < b.size(); ++i){
    //    b[i]/= norm_b;
}

void run_unconstrained_error_with(int N, int d, int nExp,double eps){
    int dim_b = nMonomials(N,2*d);
    vector<double> b(dim_b,0);
    randomPolynomialSolution_SOS_Solution0(b,N,d);
    vector<vector<int>> alpha = generate_all_exponents(N,2*d);

    double gamma_0 = 3;
    PolyTable polytable;
    get<0>(polytable) = N;
    get<1>(polytable) = 0;
    get<2>(polytable) = 0;
    get<7>(polytable) = alpha;
    get<6>(polytable) = b;
    sdp SDP(polytable,d);
    admm ADMM(SDP);
    tuple<double,int,string,timestamp_t,double,double> sol;
    vector<double> objective(nExp);
    vector<int> iterations(nExp);
    vector<timestamp_t> times(nExp);
    vector<double> pres(nExp);
    vector<double> dres(nExp);
    vector<int> errorCases(7); // < 1e-5 , [1e-5,1e-4],[1e-4,1e-3],[1e-3,1e-2],[1e-2,1e-1],[1e-1,1],>1
    vector<double> timesCases(10); // Cutoff 1000s instead of 1
    vector<double> primRes(nExp);
    vector<double> dualRes(nExp);
    for(int exp = 0; exp < nExp; ++exp){
        randomPolynomialSolution_SOS_Solution0(b,N,d);
        for(int i = 0; i < b.size(); ++i){
            b[i] *= 2;
        }
        b[0] = gamma_0;
        get<6>(polytable) = b;
        sol = callAdmm(polytable,d,eps);
        objective[exp] = get<0>(sol)-gamma_0;
        iterations[exp] = get<1>(sol);
        times[exp] = get<3>(sol);
        pres[exp] = get<4>(sol);
        dres[exp] = get<5>(sol);
    }
    // Evaluate
    for(int i = 0; i < nExp; ++i){
        if(log10(abs(objective[i])) < -5) ++errorCases[0];
        else if(log10(abs(objective[i])) > 0) ++errorCases[6];
        else ++errorCases[int(log10(abs(objective[i])))+5];
    }
    double time_i;
    for(int i = 0; i < nExp; ++i){
        time_i = double(times[i])/1e6;
        if(log10(time_i) < -5) ++timesCases[0];
        else if(log10(time_i) > 3) ++timesCases[9];
        else ++timesCases[int(log10(time_i))+5];
    }
    cout << "N=" << N << " d=" << d << " : ";
    abs(objective);
    cout << "Err=" << avg(objective) << " | ";
    cout << "Time=" << avg(times)/1e3 << " ms | ";
    cout << "Iter=" << avg(iterations) << " | ";
    printVectorNoNewLine(timesCases,"TimeCases");
    cout << "|";
    printVectorNoNewLine(errorCases,"ErrorCases");
    cout << endl;
    //printVector(pres,"pres");
    //printVector(dres,"dres");

}

void  run_unconstrained_error() {
    //vector<double> u(5,2);
    //vector<double> x(5,7);
    //cout << normOfDiffsq(u,x) << endl;
    //cout << maximum(1,2,3,4) << endl;
    //cout << norm(u) << endl;
    //pres = sqrt(normOfDiffsq(u,x)+normOfDiffsq(HiX,z));
    //double pnorm = maximum(norm(u),norm(x),norm(HiX),norm(z));
    //return;
    for (int N = 2; N <= 6; N += 2) {
        for (double eps = 1e-4; eps >= 1e-10; eps /= 10) {
            for (int d = 2; d <= 2; d += 1) {
                run_unconstrained_error_with(N, d, 50, eps);
            }
        }
    }
}



#endif //SOSADMM_BOLT_UNCONSTRAINEDTEST_H
