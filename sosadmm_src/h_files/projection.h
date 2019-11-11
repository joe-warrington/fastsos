//
// Created by Jan Schappi on 29.10.19.
//

#ifndef UNTITLED5_LAPACKPROJECTION_H
#define UNTITLED5_LAPACKPROJECTION_H

#include <vector>
#include <iostream>
#include "../../timing.h"
using namespace std;

//void projectOnPSD_LAPACK(vector<double>& x, int N);
void useLAPACK(vector<double>& x, int N);
void projectOnCone_LAPACK(vector<double>& x, int col, int coneSize);
void useLAPACK(vector<double>& x, int offset, int N);
void VLVT(vector<double>& x, double* w, double* z, int m, int N);
void writeBack(vector<double>& x, double* z, double* w, int m, int N);

#endif //UNTITLED5_LAPACKPROJECTION_H