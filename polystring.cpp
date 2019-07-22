// (c) 2019 ETH Zurich, Automatic Control Lab, Joe Warrington

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <string>

using namespace std;

struct PolyInfo{
    int dimension = 0;
    int degree = 0;
    int n_terms = 0;
};

PolyInfo infer_dim_and_n_terms(const string& s, bool diag_msg) {
    PolyInfo output;

    int n = 1; // Dimension (even constants should be assumed to live in a one-dimensional space of variables)
    int d = 0; // Degree
    int t = 0; // Number of terms found

    string number_string;
    bool still_digits;
    int inc = 0;

    bool expect_coefficient = true;
    bool expect_subscript = false;
    bool expect_exponent = false;
    bool found_exponent = false;

    int current_dimension = 0;
    int degree_this_monomial = 0;

    if (diag_msg) cout << "  " << s << endl;

    for (string::size_type j = 0; j < s.size(); ++j) {
//        if (s[j] == ' ' || s[j] == '\n') {
//            continue;  // Completely ignore spaces and carriage returns
//        } else
        if (s[j] == 'x') {
            if (j == 0) {t += 1;}
            // New x
            if (!found_exponent && current_dimension != 0)
                degree_this_monomial += 1;
            expect_coefficient = false;
            expect_subscript = true;
            expect_exponent = false;
            current_dimension = 0;
        } else if (s[j] == '^') {
            // Exponent
            expect_exponent = true;
            expect_subscript = false;
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
                found_exponent = true;
                expect_exponent = false;
            } else if (expect_subscript) {
                if (diag_msg) {cout << "Found subscript: " << stoi(number_string) << endl;}
                n = max(n, stoi(number_string));
                current_dimension = stoi(number_string);
                expect_subscript = false;
                found_exponent = false;
            }
            j += inc - 1;  // Skip ahead inc positions in overall string
        } else if (s[j] == '+' || s[j] == '-') {
            // A plus or minus indicates a new term (i.e. start of a new monomial)
            if (!found_exponent && current_dimension != 0)
                degree_this_monomial += 1;
            t += 1;
            d = max(d, degree_this_monomial);
            expect_coefficient = true;
            expect_subscript = false;
            expect_exponent = false;
            found_exponent = false;
            degree_this_monomial = 0;
            current_dimension = 0;
        }
    }
    if (!found_exponent && current_dimension != 0)
        degree_this_monomial += 1;
    d = max(d, degree_this_monomial);
    output.dimension = n;
    output.degree = d;
    output.n_terms = t;

    if (diag_msg) {
        cout << "Dimension " << output.dimension << ", degree: " << output.degree << ", " << output.n_terms << " terms." << endl;
    }

    return output;
}

void print_polynomial_table(vector<double>& mono_coeffs, vector<vector<int> >& mono_exponents) {
    int n_terms = mono_coeffs.size();
    int dimension = mono_exponents[0].size();
    for (int t = -1; t < n_terms; ++t) {
        if (t == -1) {
            cout << "coeff \t|";
            for (int i = 0; i < dimension; ++i) {
                cout << "x" << i + 1 << "\t";
            }
        } else {
            cout << mono_coeffs[t] << "\t|";
            for (int i = 0; i < dimension; ++i) {
                cout << mono_exponents[t][i] << "\t";
            }
        }
        cout << endl;
    }
}

void parse_poly(string& s, vector<double>& mono_coeffs, vector<vector<int> >& mono_exponents, PolyInfo& info,
        bool diag_msg) {

    int t = 0; // Term counter

    string number_string;
    bool still_digits;
    int inc = 0;

    int coeff_sign = 1;
    bool expect_coefficient = true;
    bool expect_x = false;
    bool expect_subscript = false;
    bool expect_exponent = false;

    bool found_coefficient = false;
    bool found_subscript = false;
    bool found_exponent = false;

    int current_dimension = 0;

    for (string::size_type j = 0; j < s.size(); ++j) {
//        if (s[j] == ' ' || s[j] == '\n') {
//            continue;  // Completely ignore spaces and carriage returns
//        } else
        if (j == 0 && s[j] != '+' && s[j] != '-')
            t = 1;  // If first character is a + or -, term will be incremented in the relevant case below
        if (s[j] == 'x') {
            if (expect_subscript) {
                cout << "Did not find a subscript for this x variable!" << endl;
                return;
            }
            // New x
            if (!found_coefficient && expect_coefficient) {
                mono_coeffs[t - 1] = coeff_sign * 1.0;  // Coefficient is 1, with sign governed by last sign char seen
            }
            if (!found_exponent && current_dimension != 0) {
                mono_exponents[t - 1][current_dimension - 1] += 1;
            }
            expect_coefficient = false;
            expect_subscript = true;
            found_subscript = false;
            expect_exponent = false;
            found_exponent = false;
            expect_x = false;
        } else if (s[j] == '^') {
            if (expect_subscript) {
                cout << "Did not find a subscript for this x variable!" << endl;
                return;
            }
            // Exponent
            expect_exponent = true;
            found_exponent = false;
            expect_subscript = false;
            found_subscript = false;
        } else if (isdigit(s[j]) || s[j] == '.') {
            // Digit is part of a coefficient, subscript or a superscript
            number_string = s[j];
            still_digits = true;
            inc = 1;
            // Collect number string
            while (still_digits) {
                if (j < s.size() - inc && (isdigit(s[j + inc]) || s[j + inc] == '.')) {
                    number_string += s[j + inc];
                    inc++;
                } else {still_digits = false;}
            }
            // Handle number string as coefficient, subscript, or exponent
            if (expect_coefficient) {
                if (diag_msg) {cout << "Found coefficient: " << stof(number_string) << endl;
                cout << "Setting coefficient for term " << t << " to " << coeff_sign * stof(number_string) << endl;}
                mono_coeffs[t - 1] = coeff_sign * stof(number_string);
                found_coefficient = true;
                expect_coefficient = false;
                expect_x = true;
            } else if (expect_exponent) {
                if (diag_msg) {cout << "Found exponent: " << stoi(number_string) << endl;
                cout << "Adding " << stoi(number_string) << " to term " << t << ", column " << current_dimension - 1 << endl;}
                if (current_dimension != 0)
                    mono_exponents[t - 1][current_dimension - 1] += stoi(number_string);
                found_exponent = true;
                expect_exponent = false;
            } else if (expect_subscript) {
                if (diag_msg) {cout << "Found subscript: " << stoi(number_string) << endl;}
                current_dimension = stoi(number_string);
                found_subscript = true;
                expect_subscript = false;
            }
            j += inc - 1;  // Skip ahead inc positions in overall string
        } else if (s[j] == '+' || s[j] == '-') {
            if (expect_subscript) {
                cout << "Did not find a subscript for this x variable!" << endl;
                return;
            }
            // A plus or minus indicates a new term (i.e. start of a new monomial)
            if (!found_exponent && current_dimension != 0)
                mono_exponents[t - 1][current_dimension - 1] += 1;
            t += 1;
            if (s[j] == '-') {coeff_sign = -1;} else {coeff_sign = 1;}
            expect_x = false;
            expect_coefficient = true;
            found_coefficient = false;
            current_dimension = 0;
        }
    }
    // Close-off actions at end of string
    if (found_coefficient && expect_x) found_exponent = true;  // Must have only found a constant at the end, so we know exponent = 0
    if (!found_exponent && found_subscript && current_dimension != 0) {
        mono_exponents[t - 1][current_dimension - 1] += 1;
    }

    if (diag_msg) {
        print_polynomial_table(mono_coeffs, mono_exponents);
    }
}

tuple<string, vector<string>, vector<string>, bool> read_problem_from_file(const string& filename) {
    string line, obj_string;
    vector<string> ineq_constr_strings(0, "");
    vector<string> eq_constr_strings(0, "");
    bool success = false;

    int f_count = 0, g_count = 0, h_count = 0;

    ifstream inFile;
    inFile.open(filename);
    if (inFile.is_open())
    {
        int row_number = 0;
        while ( getline (inFile, line) ) {
            if (line.substr(0, 2).compare("f ") == 0) {
                obj_string = line.substr(2);  // Read objective function from first line, characters 2 to end
                f_count++;
            } else if (line.substr(0, 2).compare("g ") == 0) {
                ineq_constr_strings.push_back(
                        line.substr(2));  // Read constraint function into vector of constraint strings
                g_count++;
            } else if (line.substr(0, 2).compare("h ") == 0) {
                eq_constr_strings.push_back(
                        line.substr(2));  // Read constraint function into vector of constraint strings
                h_count++;
            } else {
                cout << "Ignored line " << row_number << " of file " << filename << ":\n\t" << line << endl;
            }
            row_number++;
        }
        inFile.close();
    }

    if (f_count == 1 && g_count == ineq_constr_strings.size() && h_count == eq_constr_strings.size()) {
        success = true;
        cout << "Problem read from " << filename << " has "
             << g_count << " inequality constraints, and " << h_count << " equality constraints.\n";
        }
    else
        cout << "WARNING: input file " << filename << " contained " << f_count << " objective terms, "
        << g_count << " inequality constraints, and " << h_count << " equality constraints.\n";

    return make_tuple(obj_string, ineq_constr_strings, eq_constr_strings, success);
}