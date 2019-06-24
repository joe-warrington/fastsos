//
// Created by Joe Warrington on 2019-05-31.
//

#include <iostream>
#include <vector>

//#include "common_structs.h"

using namespace std;

struct PolyInfo{
    int dimension = 0;
    int degree = 0;
    int n_terms = 0;
};

PolyInfo infer_dim_and_n_terms(const string& s) {
    PolyInfo output;

    int n = 0; // Dimension
    int d = 0; // Degree
    int t = 0; // Number of terms found

    string number_string;
    bool still_digits;
    int inc = 0;

    bool expect_coefficient = true;
    bool expect_subscript = false;
    bool expect_exponent = false;

    int degree_this_monomial = 0;

    bool diag_msg = false;

    for (string::size_type j = 0; j < s.size(); ++j) {
//        if (s[j] == ' ' || s[j] == '\n') {
//            continue;  // Completely ignore spaces and carriage returns
//        } else
        if (s[j] == 'x') {
            if (j == 0) {t += 1;}
            // New x
            expect_coefficient = false;
        } else if (s[j] == '^') {
            // Exponent
            expect_exponent = true;
            expect_subscript = false;
        } else if (s[j] == '_') {
            // Subscript
            expect_subscript = true;
            expect_exponent = false;
        } else if (isdigit(s[j]) || s[j] == '.') {
            if (j == 0) {t += 1;}
            // Digit is part of a coefficient, subscript or a superscript
            number_string = s[j];
            still_digits = true;
            inc = 1;
            while (still_digits) {
                if (j < s.size() - inc && (isdigit(s[j + inc]) || s[j + inc] == '.')) {
                    number_string += s[j + inc];
                    inc++;
                } else {still_digits = false;}
            }
            if (expect_coefficient){
                if (diag_msg) {cout << "Found coefficient: " << stof(number_string) << endl;}
                expect_coefficient = false;
            } else if (expect_exponent) {
                if (diag_msg) {cout << "Found exponent: " << stoi(number_string) << endl;}
                degree_this_monomial += stoi(number_string);
                expect_exponent = false;
            } else if (expect_subscript) {
                if (diag_msg) {cout << "Found subscript: " << stoi(number_string) << endl;}
                n = max(n, stoi(number_string));
                expect_subscript = false;
            }
            j += inc - 1;  // Skip ahead inc positions in overall string
        } else if (s[j] == '+' || s[j] == '-') {
            // A plus or minus indicates a new term (i.e. start of a new monomial)
            t += 1;
            d = max(d, degree_this_monomial);
            degree_this_monomial = 0;
            expect_coefficient = true;
        }
    }
    d = max(d, degree_this_monomial);
    output.dimension = n;
    output.degree = d;
    output.n_terms = t;

    if (diag_msg) {
        cout << "Dimension " << output.dimension;
        cout << ", degree: " << output.degree;
        cout << ", no. of terms: " << output.n_terms << endl;
    }

    return output;
}

void parse_poly(string& s, vector<double>& mono_coeffs, vector<vector<int> >& mono_exponents, PolyInfo& info) {

    int t = 0; // Term counter

    bool print_diag = false;

    string number_string;
    bool still_digits;
    int inc = 0;

    bool expect_coefficient = true;
    int coeff_sign = 1;
    bool expect_subscript = false;
    bool expect_exponent = false;
    bool found_exponent = true;
    bool found_coefficient = false;
    int current_dimension = 1;

    bool diag_msg = false;

    for (string::size_type j = 0; j < s.size(); ++j) {
//        if (s[j] == ' ' || s[j] == '\n') {
//            continue;  // Completely ignore spaces and carriage returns
//        } else
        if (s[j] == 'x') {
            if (j == 0) {t += 1;}
            // New x
            if (!found_coefficient) {
                mono_coeffs[t - 1] = coeff_sign * 1;
            }
            expect_coefficient = false;
            if (!found_exponent) {
                mono_exponents[t - 1][current_dimension - 1] = 1;
            }
            found_exponent = false;
        } else if (s[j] == '^') {
            // Exponent
            expect_exponent = true;
            expect_subscript = false;
        } else if (s[j] == '_') {
            // Subscript
            expect_subscript = true;
            expect_exponent = false;
        } else if (isdigit(s[j]) || s[j] == '.') {
            if (j == 0) {t += 1;}
            // Digit is part of a coefficient, subscript or a superscript
            number_string = s[j];
            still_digits = true;
            inc = 1;
            while (still_digits) {
                if (j < s.size() - inc && (isdigit(s[j + inc]) || s[j + inc] == '.')) {
                    number_string += s[j + inc];
                    inc++;
                } else {still_digits = false;}
            }
            if (expect_coefficient){
                if (diag_msg) {cout << "Found coefficient: " << stof(number_string) << endl;
                cout << "Setting coefficient for term " << t << " to " << coeff_sign * stof(number_string) << endl;}
                mono_coeffs[t - 1] = coeff_sign * stof(number_string);
                found_coefficient = true;
                expect_coefficient = false;
            } else if (expect_exponent) {
                if (diag_msg) {cout << "Found exponent: " << stoi(number_string) << endl;
                cout << "Adding " << stoi(number_string) << " to term " << t << ", column " << current_dimension - 1 << endl;}
                mono_exponents[t - 1][current_dimension - 1] = stoi(number_string);
                found_exponent = true;
                expect_exponent = false;
            } else if (expect_subscript) {
                if (diag_msg) {cout << "Found subscript: " << stoi(number_string) << endl;}
                current_dimension = stoi(number_string);
                expect_subscript = false;
            }
            j += inc - 1;  // Skip ahead inc positions in overall string
        } else if (s[j] == '+' || s[j] == '-') {
            // A plus or minus indicates a new term (i.e. start of a new monomial)
            if (!found_exponent) {mono_exponents[t - 1][current_dimension - 1] = 1;}
            t += 1;
            expect_coefficient = true;
            found_coefficient = false;
            if (s[j] == '-') {coeff_sign = -1;} else {coeff_sign = 1;}
        }
    }
    if (!found_exponent) {mono_exponents[t - 1][current_dimension - 1] = 1;}

    if (print_diag) {
        for (t = -1; t < info.n_terms; ++t) {
            if (t == -1) {
                cout << "coeff \t|";
                for (int i = 0; i < info.dimension; ++i) {
                    cout << "x" << i + 1 << "\t";
                }
            } else {
                cout << mono_coeffs[t] << "\t|";
                for (int i = 0; i < info.dimension; ++i) {
                    cout << mono_exponents[t][i] << "\t";
                }
            }
            cout << endl;
        }
    }
}
