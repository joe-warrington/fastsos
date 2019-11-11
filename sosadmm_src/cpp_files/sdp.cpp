//
// Created by Jan Schappi on 02.11.19.
//

#include "sdp.h"

sdp::sdp(PolyTable polytable, int d) {

    N = get<0>(polytable);
    m = nMonomials(N,2*d);
    this->d = d;

    // Count the number of entries
    vector<int> filledA(m,0);
    n = populateA(polytable,filledA,true);

    rowPointA = vector<int>(m+1);
    fillRowPoint(filledA);
    NNZ = rowPointA[m];
    colsA = vector<int>(NNZ);
    valsA = vector<double>(NNZ);

    populateA(polytable,filledA,false);

    b = vector<double>(m,0);
    populateB(polytable);

    colsC = vector<int>(1);
    valsC = vector<double>(1);
    populateC();

    Dinv = vector<double>(n,0);
    populateDinv();

    HiaiSqInv = vector<double>(m,0);
    populateHiaiSqInv();
}

void sdp::populateHiaiSqInv() {

    int NNZ_row;
    double* pos = &valsA[0];
    for(int row = 0; row < m; ++row){
        NNZ_row = NNZi(row);
        for(int j = 0; j < NNZ_row; ++j){
            HiaiSqInv[row] += *(pos)* *(pos);
            ++pos;
        }
    }
    for(int i = 0; i < HiaiSqInv.size(); ++i){
        HiaiSqInv[i] = 1/HiaiSqInv[i];
    }
}

void sdp::populateDinv() {

    for(int i = 0; i < colsA.size(); ++i){
        ++Dinv[colsA[i]];
    }
    for(int j = 0; j < Dinv.size(); ++j){
        Dinv[j] = 1/(1+Dinv[j]);
    }
}

void sdp::subC(vector<double> &x, double rho) {

    for(int i = 0; i < colsC.size(); ++i){
        x[colsC[i]] -= valsC[i]/rho;
    }
}

void sdp::populateC(){

    colsC[0] = 0;
    valsC[0] = -1;

}

void sdp::populateB(PolyTable& polytable){

    int nObjMon = get<6>(polytable).size();
    int row;
    for(int iMon = 0; iMon < nObjMon; ++iMon){
        row = positionOfAlpha(get<7>(polytable)[iMon]);
        b[row] = get<6>(polytable)[iMon];
    }

}

int sdp::populateA(PolyTable& polytable, vector<int>& filledA, bool count){

    int col = 0;
    col = gammaSlice(filledA,col,count);
    col = relaxationSlice(filledA,col,count);

    int nIneq = get<1>(polytable);
    int di;
    for(int iIneq = 0; iIneq < nIneq; ++iIneq){
        PolyInfo info = get<4>(polytable)[iIneq];
        di = (2*d-info.degree)/2;
        col = inequalitySlice(filledA,col,get<8>(polytable)[iIneq],get<9>(polytable)[iIneq],N,di,count);
    }

    int nEq = get<2>(polytable);
    for(int iEq = 0; iEq < nEq; ++iEq){
        PolyInfo info = get<5>(polytable)[iEq];
        di = 2*d-info.degree;
        col = equalitySlice(filledA,col,get<10>(polytable)[iEq],get<11>(polytable)[iEq],N,di, count);
    }
    return col;
}


void sdp::push_back_A(vector<int>& filledA, double val, int row, int col, bool count) {
    if(!count){
        int pos = rowPointA[row]+filledA[row];
        valsA[pos] = val;
        colsA[pos] = col;
    }
    ++filledA[row];
}

int sdp::gammaSlice(vector<int>& filledA, int col, bool count) {
    push_back_A(filledA,1, 0, col,count);
    if(count) NPSD_List.push_back({col,col+1});
    return col+1;
}

int sdp::relaxationSlice(vector<int> &filledA, int col, bool count) {

    int colBegin = col;
    // Travel imaginary Q-Matrix and add to A (offset 1)
    vector<vector<int>> alphas = generate_all_exponents(N, d);
    vector<int> alpha_lk(N);
    int row;
    int coneSize = nMonomials(N,d);

    for(int k=0; k<coneSize; ++k){
        for(int l=0; l<coneSize; ++l){
            addSmall(alpha_lk,alphas[l],alphas[k]); //At (l,k) find the alpha belonging to the entry
            row = positionOfAlpha(alpha_lk); // Map alpha_lk to i in polynomials of degree 2d (!!)
            push_back_A(filledA,1,row,col,count);
            ++col;
        }
    }
    if(count) PSD_List.push_back({colBegin,col,coneSize});
    return col;
}

int sdp::inequalitySlice(vector<int> &filledA, int col, vector<double>& polyCoeffs, vector<vector<int>>& polyAlphas, int N, int di, bool count) {

    int colBegin = col;
    vector<int> repDimensions = representedDimensions(polyAlphas);
    vector<vector<int>> padAlpha = paddedAlpha(repDimensions,N,di);
    int coneSize = padAlpha.size();

    vector<int> alpha0(N), alphaRow(N);
    int row;
    for(int k = 0; k < padAlpha.size(); ++k){
        for(int l = 0; l < padAlpha.size(); ++l) {
            addSmall(alpha0,padAlpha[k],padAlpha[l]);
            for(int polyIt = 0; polyIt < polyCoeffs.size(); ++polyIt){
                addSmall(alphaRow,polyAlphas[polyIt],alpha0);
                row = positionOfAlpha(alphaRow);
                push_back_A(filledA,polyCoeffs[polyIt],row,col,count);
            }
            ++col;
        }
    }

    if(count) PSD_List.push_back({colBegin,col,coneSize});
    return col;
}

int sdp::equalitySlice(vector<int> &filledA, int col, vector<double>& polyCoeffs, vector<vector<int>>& polyAlphas, int N, int di, bool count){

    int colBegin = col;
    vector<int> repDimensions = representedDimensions(polyAlphas);
    vector<vector<int>> padAlpha = paddedAlpha(repDimensions,N,di);

    int row;
    vector<int> alphaRow(N);
    for(int l = 0; l < padAlpha.size(); ++l){
        for(int polyIt = 0; polyIt < polyCoeffs.size(); ++polyIt){
            addSmall(alphaRow,padAlpha[l],polyAlphas[polyIt]);
            row = positionOfAlpha(alphaRow);
            push_back_A(filledA,polyCoeffs[polyIt],row,col,count);
        }
        ++col;
    }
    if(count) NPSD_List.push_back({colBegin,col});
    return col;
}


vector<vector<int>> sdp::paddedAlpha(vector<int>& representedDimensions, int N, int di){
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

vector<int> sdp::representedDimensions(const vector<vector<int>>& polyAlphas){
    vector<int> representList;
    int N = polyAlphas[0].size();
    for(int i = 0; i < N; ++i){
        for(int i_alpha = 0; i_alpha < polyAlphas.size(); ++i_alpha){
            if(polyAlphas[i_alpha][i] != 0){
                representList.push_back(i);
                break;
            }
        }
    }
    return representList;
}

void sdp::fillRowPoint(vector<int>& filledA) {
    rowPointA[0]=0;
    for(int i = 1; i < rowPointA.size(); ++i){
        rowPointA[i] += rowPointA[i-1]+filledA[i-1];
        filledA[i-1] = 0;
    }
}

const int sdp::NNZi(int row) {
    return rowPointA[row+1]-rowPointA[row];
}

void sdp::print(){
    cout << "------------------------" << endl;
    cout << "Printing SDP" << endl << endl;

    printVector(valsA,"valsA");
    printVector(colsA,"colsA");
    printVector(rowPointA,"rowPointA");

    int pos;
    cout << endl << "c = ";
    pos = 0;
    for(int l = 0; l < colsC.size(); ++l){
        for(int k = pos; k < colsC[l]; ++k){
            cout << "0 ";
        }
        pos = colsC[l]+1;
        cout << valsC[l] << " ";
    }
    for(int l = pos; l < n; ++l) cout << "0 ";
    cout << endl << endl;

    for(int i = 0; i < m; ++i){
        if(i == m/2) cout << "A = ";
        else cout << "    ";
        pos = 0;
        for(int l = 0; l < NNZi(i); ++l){
            for(int k = pos; k < colsA[rowPointA[i]+l]; ++k){
                cout << "0 ";
            }
            pos = colsA[rowPointA[i]+l]+1;
            cout << valsA[rowPointA[i]+l] << " ";
        }
        for(int l = pos; l < n; ++l) cout << "0 ";

        if(i == m/2) cout << "  b = ";
        else cout << "      ";
        cout << b[i] << endl;
    }
    cout << endl;
    cout << "Dinv = ";
    for(int i = 0; i < n; ++i){
        cout << Dinv[i] << " ";
    }
    cout << endl << "HiaiSqInv = ";
    for(int i = 0; i < m; ++i){
        cout << HiaiSqInv[i] << " ";
    }

    cout << endl << endl;
    cout << "PSD_List; {first,last+1,coneSize}" << endl;
    for(int iPSD = 0; iPSD < PSD_List.size(); ++iPSD){
        cout << iPSD << "|" << PSD_List[iPSD][0] << "," << PSD_List[iPSD][1] << "," << PSD_List[iPSD][2] << endl;
    }

    cout << endl;
    cout << "NPSD_List; {first,last+1}" << endl;
    for(int iNPSD = 0; iNPSD < NPSD_List.size(); ++iNPSD){
        cout << iNPSD << "|" << NPSD_List[iNPSD][0] << "," << NPSD_List[iNPSD][1] << endl;
    }

    cout << endl;
    cout << "m=" << m;
    cout << "  n=" << n;
    cout << "  NNZ=" << NNZ;
    cout << "  d=" << d;
    cout << "  N=" << N << endl;
}