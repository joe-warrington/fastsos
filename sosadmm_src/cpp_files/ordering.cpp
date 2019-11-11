//
// Created by Jan Schappi on 02.11.19.
//

#include "ordering.h"

int positionOfAlpha(const vector<int>& alpha){
    int n = alpha.size();
    int q; int i = 0; int rel_i;
    for(int start = 0; start < n; ++start){
        q = sum(alpha,start,n);
        if(q == 0) return i;
        else{
            rel_i = n_choose_k(n-start+q-1,n-start);
            i += rel_i;
        }
    }
    return i;
}

int sum(const vector<int>& alpha,const int start,const int end){
    int s = 0;
    for(int i = start; i < alpha.size() && i < end; ++i){
        s += alpha[i];
    }
    return s;
}

unsigned long int n_choose_k(int a, int b){
    // Formula for a-choose-b is a!/[b!(a-b)!
    if(a <= 0 || b < 0) {
        //cout << COLOR_RED << "Cannot compute a-choose-b for a = " << a << " and b = " << b << "!" << endl;
        exit(1);
    }
    else if (b > a) {
        //cout << COLOR_RED << "Cannot compute a-choose-b for b greater than a! (" << b << " > " << a << ")" << endl;
        exit(1);
    }
    else if (b == a) {
        return 1;
    }
    else {  // Must have a > b > 0
        unsigned long int num = 1, den = 1;
        for (int i = b + 1; i <= a; i++) {
            num *= i;
        } // generates a! / b! for a > b, skipping all the cancelled terms 1 * 2 * ... * b
        for (int i = 1; i <= a - b; i++) {
            den *= i;
        } // generates (a - b)! in the denominator of the a-choose-b formula
        return num / den;
    }
}

unsigned long int nMonomials(int N, int d){
    unsigned long int result;
    int choosing_choice = max(N,d);
    // n+d choose d is the same as n+d choose n, but better overflow characteristics when the larger is chosen
    result = n_choose_k(N + d, choosing_choice);
    return result;
}

vector<vector<int>> generate_all_exponents(int N, int d){
    unsigned long int s_of_d = nMonomials(N, d);
    vector<int> blank_row(N, 0);
    auto row = blank_row;
    vector<vector <int> > vec_out(s_of_d, blank_row);
    int current_d = 1;
    vector<int> e_posns(d, 0);  // Positions of the up to d exponents
    bool new_d = true;  // Whether d was just incremented
    for (unsigned long int i = 1; i < s_of_d; i++) {
        // i = 0 corresponds to d = 0 and a row of all zeros, as initialized. Can therefore start from i = 1 (d = 1).
        row = blank_row;
        if (new_d) {
            for (int j = 0; j < current_d; j++)
                e_posns[j] = 0;  // Set position of all the current_d "ones" to 0.
            // Only the first (current_d) elements of e_posns are used.
            new_d = false;
        }
        else {
            if (e_posns[current_d - 1] == N - 1 && current_d > 1) {
                for (int j = current_d - 2; j >= 0; j--) {
                    // Search back through more significant "ones" to find one that hasn't reached the end posn (n - 1)
                    if (e_posns[j] < N - 1) {
                        e_posns[j] += 1;
                        for (int k = j + 1; k < current_d; k++)
                            e_posns[k] = e_posns[j];  // Set positions of all less significant "ones" to that of jth.
                        break;
                    }
                }
            } else {
                e_posns[current_d - 1] += 1;  // Advance least significant "one"
            }
        }
        for (int j = 0; j < current_d; j++)
            row[e_posns[j]] += 1;  // Add the exponents to the correct columns of the row
        vec_out[i] = row;  // Add the row to the output.

        if (row[N - 1] == current_d && i != s_of_d - 1) {
            new_d = true;
            current_d += 1;
        }
    }
    return vec_out;
}
