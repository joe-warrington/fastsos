//
// Created by Jan Schappi on 08.11.19.
//

#ifndef SOSADMM_BOLT_CONSTRAINEDTEST_H
#define SOSADMM_BOLT_CONSTRAINEDTEST_H

#include "sdp.h"
#include "unconstrainedTest.h"

void printPolynomial(vector<double> coeffs, vector<vector<int>> alphas){
    bool firstPrintDone = false;
    int N = alphas[0].size();
    for(int i = 0; i < coeffs.size(); ++i){
        if(i != 0 && coeffs[i] > 0 && firstPrintDone) cout << "+ ";
        if(coeffs[i] != 0) {
            if(coeffs[i] > 0) cout << setprecision(3) << coeffs[i];
            else cout << " - " << setprecision(3) << -coeffs[i];
            firstPrintDone = true;
            for (int j = 0; j < N; ++j) {
                if (alphas[i][j] != 0) {
                    cout << "x" << j+1;
                    if (alphas[i][j] != 1) cout << "^" << alphas[i][j];
                }
            }
            cout << " ";
        }
    }
    cout << endl;
}

// Samples k random integers in [N]
vector<int> randomDimensions(int N, int k){
    vector<int> notPicked(N);
    vector<int> picked(k);
    int randInt;
    for(int i = 0; i < N; ++i) notPicked[i] = i;
    for(int i = 0; i < k; ++i){
        randInt = rand() % (N-i);
        picked[i] = notPicked[randInt];
        notPicked.erase(notPicked.begin()+randInt);
    }
    return picked;
}

vector<vector<int>> paddedAlpha(vector<int>& representedDimensions, int N, int di){
    int Ni = representedDimensions.size();
    vector<vector<int>> AlphaReduced = generate_all_exponents(Ni,di);
    vector<vector<int>> AlphaPadded;
    int dim;
    for(int i = 0; i < AlphaReduced.size(); ++i){
        vector<int> alpha(N,0);
        AlphaPadded.push_back(alpha);
        for(int repDimIt = 0; repDimIt < representedDimensions.size(); ++repDimIt){
            dim = representedDimensions[repDimIt];
            AlphaPadded[i][dim] = AlphaReduced[i][repDimIt];
        }
    }
    return AlphaPadded;
}

// coeff is already initialized to 0 and has length nMonomials(N-k,d)
void sparseRandomPolynomials(vector<double>& coeff, vector<vector<int>>& alphas, int N, int d, int k){
    vector<int> toUse = randomDimensions(N,k);
    totallyRandomPolynomial(coeff,alphas,k,d);
    alphas = paddedAlpha(toUse,N,d);
    coeff.resize(alphas.size());
}

void printPolyTable(PolyTable& polytable, string name){
    cout << "----------------------" << endl;
    cout << "PolyTable " << name << endl;
    cout << "Dimensions: " << get<0>(polytable) << endl;
    cout << "Inequalities: " << get<1>(polytable) << endl;
    cout << "Equalities: " << get<2>(polytable) << endl;
    cout << "PolyInfo (Objective): " << " degree=" << get<3>(polytable).degree << " dimension=" << get<3>(polytable).dimension << " nterms=" << get<3>(polytable).n_terms << endl;
    for(int i = 0; i<get<1>(polytable); ++i){
        cout << "PolyInfo (Ineq[" << i << "]): " << " degree=" << get<4>(polytable)[i].degree << " dimension=" << get<4>(polytable)[i].dimension << " nterms=" << get<4>(polytable)[i].n_terms << endl;
        printPolynomial(get<8>(polytable)[i],get<9>(polytable)[i]);
    }
    for(int i = 0; i<get<2>(polytable); ++i){
        cout << "PolyInfo (Eq[" << i << "]): " << " degree=" << get<5>(polytable)[i].degree << " dimension=" << get<5>(polytable)[i].dimension << " nterms=" << get<5>(polytable)[i].n_terms << endl;
        printPolynomial(get<10>(polytable)[i],get<11>(polytable)[i]);
    }
}

void printOptimizationProblem(PolyTable& polytable, string name){
    cout << "----------------------" << endl;
    cout << "Problem " << name << endl;
    cout << "f ";
    printPolynomial(get<6>(polytable),get<7>(polytable));
    for(int i = 0; i<get<1>(polytable); ++i){
        cout << "g ";
        printPolynomial(get<8>(polytable)[i],get<9>(polytable)[i]);
    }
    for(int i = 0; i<get<2>(polytable); ++i){
        cout << "h ";
        printPolynomial(get<10>(polytable)[i],get<11>(polytable)[i]);
    }
}

tuple<double,int,int,timestamp_t> constrainedRandomSparseTest_with(int N, int dObj, int nIneq, int dIneq, int nEq, int dEq, int k){
    // Minimizing an SOS-Polynomial with solution at 0 (can be done as in unconstrained test)
    // All inequality constraints have a y-intercept > 0. Hence they are feasible at 0
    // All equality constraints have a y-intercept of 0.
    if(dObj % 2 != 0) {
        cout << "Warning: Objective degree should be even for this setup (sampling SOS)" << endl;
        exit(2);
    }
    int mObj = nMonomials(N,dObj);
    int mIneq = nMonomials(N,dIneq);
    int mEq = nMonomials(N,dEq);
    vector<double> p_coeffs(mObj);
    vector<vector<int>> p_alphas;
    vector<vector<double>> g_coeffs;
    vector<vector<vector<int>>> g_alphas;
    vector<PolyInfo> g_info;
    vector<vector<double>> h_coeffs;
    vector<vector<vector<int>>> h_alphas;
    vector<PolyInfo> h_info;

    // Objective
    randomPolynomialSolution_SOS_Solution0(p_coeffs, N, dObj/2); // Returns a random polynomial of degree 2d
    p_alphas = generate_all_exponents(N,dObj);
    PolyInfo p_info;
    p_info.degree = dObj; p_info.dimension = N; p_info.n_terms = nMonomials(N,dObj);
    for(int i = 0; i < nIneq; ++i){

        g_coeffs.push_back(vector<double>(mIneq));
        g_alphas.push_back(vector<vector<int>>(mIneq));
        sparseRandomPolynomials(g_coeffs[i],g_alphas[i],N,dIneq,k);
        g_coeffs[i][0] = 1;
        PolyInfo info; info.degree = dIneq; info.dimension = N; info.n_terms = g_coeffs[i].size();
        g_info.push_back(info);
    }
    for(int i = 0; i < nEq; ++i){
        h_coeffs.push_back(vector<double>(mEq));
        h_alphas.push_back(vector<vector<int>>(mEq));
        sparseRandomPolynomials(h_coeffs[i],h_alphas[i],N,dEq,k);
        h_coeffs[i][0] = 0;
        PolyInfo info; info.degree = dEq; info.dimension = N; info.n_terms = h_coeffs[i].size();
        h_info.push_back(info);
    }
    PolyTable polytable;
    get<0>(polytable)= N;
    get<1>(polytable)= nIneq;
    get<2>(polytable)= nEq;
    get<4>(polytable) = g_info;
    get<5>(polytable) = h_info;
    get<6>(polytable)= p_coeffs;
    get<7>(polytable)= p_alphas;
    get<8>(polytable)= g_coeffs;
    get<9>(polytable)= g_alphas;
    get<10>(polytable)= h_coeffs;
    get<11>(polytable)= h_alphas;
    printOptimizationProblem(polytable,"Constrained Random Problem");
    int d = (dObj+1)/2+3;
    if((dIneq+1)/2 > d) d = (dIneq+1)/2;
    if((dEq+1)/2 > d) d = (dEq+1)/2;
    cout << "Building sdp..." << endl;
    sdp SDP(polytable,d);
    cout << "sdp built. Optimizing.." << endl;
    admm ADMM(SDP);
    auto result = ADMM.min();
    cout << "Found solution " << get<0>(result) << endl;
    cout << "Iterations " << ADMM.itUsed << endl;
    cout << "Total Time " << ADMM.tTotal << endl;
    cout << get<2>(result) << endl;
    cout << "Optimization finished." << endl;
    return make_tuple(get<0>(result),get<1>(result),ADMM.itUsed,ADMM.tTotal);
}

tuple<double,int,int,timestamp_t> constrainedRandomTest_with(int N, int dObj, int nIneq, int dIneq, int nEq, int dEq){
    // Minimizing an SOS-Polynomial with solution at 0 (can be done as in unconstrained test)
    // All inequality constraints have a y-intercept > 0. Hence they are feasible at 0
    // All equality constraints have a y-intercept of 0.
    if(dObj % 2 != 0) {
        cout << "Warning: Objective degree should be even for this setup (sampling SOS)" << endl;
        exit(2);
    }
    int mObj = nMonomials(N,dObj);
    int mIneq = nMonomials(N,dIneq);
    int mEq = nMonomials(N,dEq);
    vector<double> p_coeffs(mObj);
    vector<vector<int>> p_alphas;
    vector<vector<double>> g_coeffs;
    vector<vector<vector<int>>> g_alphas;
    vector<PolyInfo> g_info;
    vector<vector<double>> h_coeffs;
    vector<vector<vector<int>>> h_alphas;
    vector<PolyInfo> h_info;

    // Objective
    randomPolynomialSolution_SOS_Solution0(p_coeffs, N, dObj/2); // Returns a random polynomial of degree 2d
    p_coeffs[0] = 3;
    p_alphas = generate_all_exponents(N,dObj);
    PolyInfo p_info;
    p_info.degree = dObj; p_info.dimension = N; p_info.n_terms = nMonomials(N,dObj);
    for(int i = 0; i < nIneq; ++i){

        g_coeffs.push_back(vector<double>(mIneq));
        g_alphas.push_back(vector<vector<int>>(mIneq));
        totallyRandomPolynomial(g_coeffs[i],g_alphas[i],N,dIneq);
        g_coeffs[i][0] = 1;
        PolyInfo info; info.degree = dIneq; info.dimension = N; info.n_terms = g_coeffs[i].size();
        g_info.push_back(info);
    }
    for(int i = 0; i < nEq; ++i){
        h_coeffs.push_back(vector<double>(mEq));
        h_alphas.push_back(vector<vector<int>>(mEq));
        totallyRandomPolynomial(h_coeffs[i],h_alphas[i],N,dEq);
        h_coeffs[i][0] = 0;
        PolyInfo info; info.degree = dEq; info.dimension = N; info.n_terms = h_coeffs[i].size();
        h_info.push_back(info);
    }
    PolyTable polytable;
    get<0>(polytable)= N;
    get<1>(polytable)= nIneq;
    get<2>(polytable)= nEq;
    get<4>(polytable) = g_info;
    get<5>(polytable) = h_info;
    get<6>(polytable)= p_coeffs;
    get<7>(polytable)= p_alphas;
    get<8>(polytable)= g_coeffs;
    get<9>(polytable)= g_alphas;
    get<10>(polytable)= h_coeffs;
    get<11>(polytable)= h_alphas;
    printOptimizationProblem(polytable,"Constrained Random Problem");
    int d = (dObj+1)/2;
    if((dIneq+1)/2 > d) d = (dIneq+1)/2;
    if((dEq+1)/2 > d) d = (dEq+1)/2;
    cout << "Building sdp..." << endl;
    sdp SDP(polytable,d);
    cout << "sdp built. Optimizing.." << endl;
    admm ADMM(SDP);
    auto result = ADMM.min();
    cout << "Found solution " << get<0>(result) << endl;
    cout << "Iterations " << ADMM.itUsed << endl;
    cout << "Total Time " << ADMM.tTotal << endl;
    cout << get<2>(result) << endl;
    cout << "Optimization finished." << endl;
    return make_tuple(get<0>(result),get<1>(result),ADMM.itUsed,ADMM.tTotal);
}

void constrainedRandomTest(){

    //printVector(coeffs,"Coeffs");
    //printVector(alphas,"Alphas");

    int N = 3;
    int dObj = 2;
    int nIneq = 0;
    int dIneq = 2;
    int nEq = 0;
    int dEq = 3;
    int nExp = 1;
    // Return
    double err_avg=0,it_avg=0,t_avg = 0;
    int validSolved = 0;
    for (int i = 0; i < nExp; ++i) {
        auto returnTuple = constrainedRandomTest_with(N, dObj, nIneq, dIneq, nEq, dEq);
        if (abs(get<0>(returnTuple)) > 1) cout << "Blunder." << endl;
        else {
            err_avg = (validSolved * err_avg + abs(get<0>(returnTuple))) / (i + 1);
            t_avg = (validSolved * t_avg + abs(double(get<3>(returnTuple)))) / (i + 1);
            it_avg = (validSolved * it_avg + abs(double(get<2>(returnTuple)))) / (i + 1);
            ++validSolved;
        }
    }
    cout << "-----------------------------------" << endl;
    cout << "Solved with error smaller than 1: " << validSolved << "/" << nExp << endl;
    cout << "Average error: " << err_avg << endl;
    cout << "Average time: " << t_avg / 1e6 << "ms" << endl;
    cout << "Average number of iterations " << it_avg << endl;

}

#endif //SOSADMM_BOLT_CONSTRAINEDTEST_H
