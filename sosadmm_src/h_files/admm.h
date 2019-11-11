//
// Created by Jan Schappi on 02.11.19.
//

#ifndef SOSADMM_BOLT_ADMM_H
#define SOSADMM_BOLT_ADMM_H

#include "sdp.h"
#include "h_ops.h"
#include "vecop.h"
#include "projection.h"
#include "../../timing.h"
#include <cmath>

struct admm{

    sdp* SDP;

    vector<double> x;
    vector<double> u;
    vector<double> w;
    vector<double> z;
    vector<double> mu;
    vector<double> zeta;
    vector<double> HiX;
    vector<double> z_temp;
    vector<double> w_temp;

    vector<double> z_old;
    vector<double> u_old;

    tuple<double,int,string> evaluation;

    double rho = 1;
    double rho_max = 1e3;
    double rho_min = 1e-2;
    double rho_mult = 1;
    int stepsPrimalSlow = 0;
    int stepsDualSlow = 0;
    int stepsOneSidedThresh = 50;

    int threads_n = 1;
    int threads_m = 1;
    int threads_NNZ = 1;

    double eps = 1e-7;
    int maxIt = 1e5;

    int checkConvergenceAfterSteps = 2;
    int checkInfeasibleAfterSteps = 1000;
    int monitorAfterSteps = 100000;
    bool monitor = true;
    bool monitorFinal = true;

    double pres = -1;
    double dres = -1;
    double presLastTime = -1;
    double dresLastTime = -1;

    timestamp_t tTotal,tIt,tP,tA,tD,tC,tPSD;

    int nThreads_n;
    int nThreads_NNZ;

    admm(sdp& SDP);
    tuple<double,int,string> min();
    void updateP();
    void updateA();
    void updateD();
    void projectOnPSD(vector<double>& u_proj);
    double obj();

    bool control(int it);
    void rhoAdaptation();
    bool convergenceCheck(int it, bool saveOld);
    int infeasibilityCheck();

    void printVariables(string status);
    void monitorStatus(int it);
    void monitorStatusFinal(int it);
    int itUsed;
};

#endif //SOSADMM_BOLT_ADMM_H
