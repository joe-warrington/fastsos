//
// Created by Jan Schappi on 29.10.19.
//

#include "projection.h"
using namespace std;

extern "C" {
extern int dsyevr_(char*,char*,char*,int*,double*,int*,double*,double*,int*,int*,
                   double*,int*,double*,double*,int*,int*,double*,int*,int*,int*,
                   int*);
}

void useLAPACK(vector<double>& x, int N){

    timestamp_t t0 = timenow();

    int* wsaved; int* zsaved; int msaved;
    {
        /* Locals */
        int n = N, il, iu, m, lda = N, ldz = N, info, lwork, liwork;
        double abstol;
        double vl, vu;
        int iwkopt;
        int *iwork;
        double wkopt;
        double *work;
        /* Local arrays */
        int isuppz[N];
        double w[N], z[N * N];

        /* Negative abstol means using the default value */
        abstol = -1.0;
        /* Set il, iu to compute NSELECT smallest eigenvalues */
        vl = 0;
        vu = 1.79769e+308;
        /* Query and allocate the optimal workspace */
        lwork = -1;
        liwork = -1;
        dsyevr_((char *) "Vectors", (char *) "V", (char *) "Upper", &n, &x[0], &lda, &vl, &vu, &il, &iu,
                &abstol, &m, w, z, &ldz, isuppz, &wkopt, &lwork, &iwkopt, &liwork,
                &info);
        lwork = (int) wkopt;
        work = (double *) malloc(lwork * sizeof(double));
        liwork = iwkopt;
        iwork = (int *) malloc(liwork * sizeof(int));
        /* Solve eigenproblem */
        dsyevr_((char *) "Vectors", (char *) "V", (char *) "Upper", &n, &x[0], &lda, &vl, &vu, &il, &iu,
                &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork,
                &info);
        /* Check for convergence */
        if (info > 0) {
            printf("The dsyevr (useLAPACK) failed to compute eigenvalues.\n");
            exit(1);
        }
        t0 = timenow() - t0;
        free((void *) iwork);
        free((void *) work);

        timestamp_t t0 = timenow();
        writeBack(x,z,w,m,N);
        t0 = timenow() - t0;
    }
}

void addvvT(vector<double>& x, double* z, double* w, int m, int N, int i){
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

void makeSymFromHalf(vector<double>& x, int N){
    for(int i = 0; i < N; ++i){
        for(int j = i+1; j < N; ++j){
            x[i+N*j] = x[j + N*i];
        }
    }

}

void writeBack(vector<double>& x, double* z, double* w, int m, int N){
    memset(&x[0],0,x.size()*sizeof(double));
    for(int i = 0; i < m; ++i){
        addvvT(x, z, w, m, N, i);
    }
    makeSymFromHalf(x,N);
}



void projectOnCone_LAPACK(vector<double>& x, int col, int N){
    vector<double> xi(N*N);
    double entry;
    for(int i = 0; i < N; ++i){
        for(int j = i; j < N; ++j){
            entry = (x[col+N*i+j] + x[col+N*j+i])/2;
            xi[N*i+j] = entry;
            xi[N*j+i] = entry;
        }
    }
    useLAPACK(xi,N);

    for(int i = 0; i < xi.size(); ++i){
        x[col+i] = xi[i];
    }
}