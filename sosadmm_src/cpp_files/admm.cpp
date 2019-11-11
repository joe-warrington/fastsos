//
// Created by Jan Schappi on 02.11.19.
//

#include "admm.h"

admm::admm(sdp& SDP) {

    x = vector<double>(SDP.n,0);
    u = vector<double>(SDP.n,0);
    w = vector<double>(SDP.m,0);
    z = vector<double>(SDP.NNZ,0);
    mu = vector<double>(SDP.NNZ,0);
    zeta = vector<double>(SDP.n,0);
    HiX = vector<double>(SDP.NNZ,0);
    z_temp = vector<double>(SDP.NNZ,0);
    w_temp = vector<double>(SDP.m,0);

    this->SDP = &SDP;
}

tuple<double,int,string> admm::min(){

    tTotal = timenow();
    int convergenceInfo;
    for(int it = 0; it < maxIt; ++it){
        if(monitor) tIt = timenow();
        updateP();
        updateA();

        if(control(it)) {
            //monitorStatus(it);
            itUsed = it;
            updateD();
            updateP();
            return evaluation;
        }

        updateD();
        tIt = timenow()-tIt;

        if(monitor) monitorStatus(it);
    }
    updateP();
    control(maxIt);
    return evaluation;
};

void admm::updateP(){
    if(monitor) tP=timenow();
    // Update Primal
    addAndMulSecond(z,z,mu,1/rho,threads_NNZ);
    sumHiT(x,z,SDP->colsA);
    addAndMulSecond(u,u,zeta,1/rho,threads_n);
    addAndMulSecond(x,x,u,1,threads_n);
    SDP->subC(x,rho);
    mulPointwise(x,SDP->Dinv,x,threads_n);
    if(monitor) tP= timenow()-tP;
}
// Modified: z,u
void admm::updateA() {

    if(monitor) tA = timenow();
    addAndMulSecond(u,x,zeta,-1/rho,threads_n);
    if(monitor) tPSD = timenow();
    projectOnPSD(u);
    if(monitor) tPSD = timenow()-tPSD;

    Hix_seq(HiX,x,SDP->colsA);
    HiaiT_seq(w,HiX,SDP->rowPointA,SDP->valsA);
    HiaiT_seq(w_temp,mu,SDP->rowPointA,SDP->valsA);
    addAndMulSecond(w,w,w_temp,-1/rho,threads_m);
    sub(w,w,SDP->b,threads_m);
    mulPointwise(w,w,SDP->HiaiSqInv,threads_m);

    addAndMulSecond(z,HiX,mu,-1/rho,threads_NNZ);
    Hiai_seq(z_temp,w,SDP->rowPointA,SDP->valsA);
    sub(z,z,z_temp,threads_NNZ);

    if(monitor) tA = timenow()-tA;
}
// Modified: None
void admm::updateD() {

    if(monitor) tD = timenow();
    sub(z_temp,z,HiX,threads_NNZ);
    addAndMulSecond(mu,mu,z_temp,rho,threads_NNZ);

    sub(x,u,x,threads_n);
    addAndMulSecond(zeta,zeta,x,rho,threads_n);

    if(monitor) tD = timenow()-tD;
}
// Modified: x
void admm::projectOnPSD(vector<double> &u_proj) {

    int coneSize; int col;
    for(int iPSD = 0; iPSD < SDP->PSD_List.size(); ++iPSD){
        col = SDP->PSD_List[iPSD][0];
        coneSize = SDP->PSD_List[iPSD][2];
        //printVector(u,"In PRJECT ON PSD");
        projectOnCone_LAPACK(u,col,coneSize);

        //printVector(u,"After proj.");
    }
}

void admm::printVariables(string status){
    cout << "-------------------------" << endl;
    cout << "ADMM (" << status << ")" << endl;
    if(x.size() < 100){
        printVector(x,"x");
        printVector(u,"u");
        printVector(w,"w");
        printVector(z,"z");
        printVector(mu,"mu");
        printVector(zeta,"zeta");
        printVector(HiX,"HiX");
        printVector(z_temp,"z_temp");
        printVector(w_temp,"w_temp");
    }

}

double admm::obj(){
    double objective = 0;
    for(int i = 0; i < SDP->valsC.size(); ++i){
        objective += (SDP->valsC[i])*(x[SDP->colsC[i]]);
    }
    return -objective;
}
bool admm::convergenceCheck(int it, bool saveOld) {

    if(saveOld){
        u_old = u;
        z_old = z;
        return false;
    }
    else{
        //pres = sqrt(normOfDiffsq(u,x)+normOfDiffsq(HiX,z));
        //double pnorm = maximum(norm(u),norm(x),norm(HiX),norm(z));
        pres = sqrt(normOfDiffsq(u,x)+normOfDiffsq(HiX,z));
        double pnorm = maximum(norm(u),norm(x),norm(HiX),norm(z));
        pres /= pnorm;

        dres = normOfDiffsq(u_old,u)+normOfDiffsq(z_old,z);
        double dnorm = normsq(mu)+normsq(zeta);
        dres = rho*sqrt(dres/dnorm);

        if(max(dres,pres) < eps) return 1;
        return 0;
    }
}
int admm::infeasibilityCheck() {

    if(presLastTime < 0){presLastTime = pres; dresLastTime = dres;}
    else{
        if(pres/dres > 1e5 && (pres-presLastTime)/pres > -0.01) return -1; // Infeasible: No progress and constr not satisfied
        if(pres/dres < 1e-5 && (dres-dresLastTime)/dres > -0.01) return -2; // Unbounded: Constraints satisfied but keeps running
    }
    return 0;
}
void admm::rhoAdaptation() {

    double rat = pres/dres;
    if(rat > 10) {
        if(stepsPrimalSlow >= stepsOneSidedThresh) rho = ::min(rho_mult*rho,rho_max);
        else{
            stepsPrimalSlow += checkConvergenceAfterSteps;
            stepsDualSlow = 0;
        }
    }
    else{
        if(stepsDualSlow >= stepsOneSidedThresh) rho = max(rho/rho_mult,rho_min);
        else{
            stepsDualSlow += checkConvergenceAfterSteps;
            stepsPrimalSlow = 0;
        }
    }
}
bool admm::control(int it){

    if(it % checkConvergenceAfterSteps == 0 && it > 0) convergenceCheck(it,true);
    if(it % checkConvergenceAfterSteps == 1 && it > 1) {
        if(monitor) tC = timenow();
        if(convergenceCheck(it,false)) {
            tTotal = timenow()-tTotal;
            double minimum = obj();
            evaluation = make_tuple(obj(),1,"Opt =" + to_string(minimum) +  "; Converged after " + to_string(it) + " iterations in " + time_string(tTotal));
            return true;
        }
        rhoAdaptation();
        if(monitor) tC = timenow()-tC;
    }
    if(it % checkInfeasibleAfterSteps == 0 && it > 0){
        int convergenceInfo = infeasibilityCheck();
        if(convergenceInfo == -2) {
            tTotal = timenow()-tTotal;
            //monitorStatusFinal(it);
            evaluation = make_tuple(obj(),convergenceInfo,"Infeasible after " + to_string(it) + " iterations in " + time_string(tTotal));
            return true;
        }
        if(convergenceInfo == -1) {
            tTotal = timenow()-tTotal;
            //monitorStatusFinal(it);
            evaluation = make_tuple(obj(),convergenceInfo,"Unbounded after " + to_string(it) + " iterations in " + time_string(tTotal));
            return true;
        }
    }
    if(it == maxIt){
        updateD();
        updateP();
        tTotal = timenow()-tTotal;
        evaluation = make_tuple(obj(),1,"Maximum number of iterations (" + to_string(maxIt) + ") exceeded. Took " + time_string(tTotal));
        return true;
    }
    return false;
}
void admm::monitorStatus(int it) {
    if(monitor && it % monitorAfterSteps == 0 && it > 0){
        cout << "---------------------" << endl;
        cout << "Monitoring (iteration " << it << ")" << endl;
        cout << "Primal Residual = " << pres << " // Dual Residual = " << dres << endl;
        cout << "Rho = " << rho << endl;
        cout << "Objective = " << obj() << endl;
        cout << "Times are " << "tP: " << time_string(tP) << " tA: " << time_string(tA) << " tD: " << time_string(tD) << " tPSD: " << time_string(tPSD) << endl;
        cout << "Percentages are " << "tP: " << 100*double(tP)/(tIt) << "% tA: " << 100*double(tA)/(tIt) << "% tD: " << 100*double(tD)/(tIt) << "% tPSD: " << 100*double(tPSD)/tIt << "%" << endl;
        cout << "Total time is " << time_string(timenow()-tTotal) << endl;
        cout << "Variable Storage: " << 4*(3*SDP->n + 4*SDP->NNZ + 2*SDP->m)/1000 << "KB" << endl;
        cout << "SDP storage: " << 4*(3*SDP->NNZ + 2*SDP->m + SDP->n)/1000 << "KB" << endl;
    }
}

void admm::monitorStatusFinal(int it) {
    cout << "---------------------" << endl;
    cout << "Monitoring Final (iteration " << it << ")" << endl;
    cout << "Primal Residual = " << pres << " // Dual Residual = " << dres << endl;
    cout << "Rho = " << rho << endl;
    cout << "Objective = " << obj() << endl;
    cout << "Times are " << "tP: " << time_string(tP) << " tA: " << time_string(tA) << " tD: " << time_string(tD) << " tPSD: " << time_string(tPSD) << endl;
    cout << "Percentages are " << "tP: " << 100*double(tP)/(tIt) << "% tA: " << 100*double(tA)/(tIt) << "% tD: " << 100*double(tD)/(tIt) << "% tPSD: " << 100*double(tPSD)/tIt << "%" << endl;
    cout << "Total time is " << time_string(timenow()-tTotal) << endl;
    cout << "Variable Storage: " << 4*(3*SDP->n + 4*SDP->NNZ + 2*SDP->m)/1000 << "KB" << endl;
    cout << "SDP storage: " << 4*(3*SDP->NNZ + 2*SDP->m + SDP->n)/1000 << "KB" << endl;
}