//
// Created by Jan Schappi on 02.11.19.
//

#ifndef SOSADMM_BOLT_SDP_H
#define SOSADMM_BOLT_SDP_H

#include <vector>
#include <array>
#include "../../manip_poly.h"
#include "ordering.h"
using namespace std;

#define PolyTable tuple<int, int, int, PolyInfo, vector<PolyInfo>, vector<PolyInfo>,vector<double>,vector<vector<int> >,vector<vector <double> >, vector<vector <vector <int> > >,vector<vector <double> >, vector<vector <vector <int> > > >

struct sdp{

    //A,b,c
    int m,n;
    vector<double> valsA;
    vector<int> colsA;
    vector<int> rowPointA;
    void push_back_A(vector<int>& filledA, double val, int row, int col, bool count);
    void fillRowPoint(vector<int>& filledA);
    const int NNZi(int row);

    vector<double> b;
    void populateB(PolyTable& polytable);

    vector<double> valsC;
    vector<int> colsC;
    void populateC();
    void subC(vector<double>& x, double rho);

    vector<double> Dinv;
    void populateDinv();
    vector<double> HiaiSqInv;
    void populateHiaiSqInv();

    //Slice Division, relaxation order and dimension (y)
    int d,N,NNZ;
    vector<array<int,3>> PSD_List; // first,last+1,coneSize
    vector<array<int,2>> NPSD_List; // first,last+1

    //Building Functions
    sdp(PolyTable polytable, int d);
    int populateA(PolyTable& polytable, vector<int>& filledA, bool count);
    int gammaSlice(vector<int>& filledA, int col, bool count);
    int relaxationSlice(vector<int>& filledA, int col, bool count);
    int inequalitySlice(vector<int> &filledA, int col, vector<double>& polyCoeffs, vector<vector<int>>& polyAlphas, int N, int di, bool count);
    int equalitySlice(vector<int> &filledA, int col, vector<double>& polyCoeffs, vector<vector<int>>& polyAlphas, int N, int di, bool count);
    vector<vector<int>> paddedAlpha(vector<int>& representedDimensions, int N, int di);
    vector<int> representedDimensions(const vector<vector<int>>& polyAlphas);

    //Print and Check
    void print();
};
#endif //SOSADMM_BOLT_SDP_H
