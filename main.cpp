#include <iostream>
#include <fstream>
#include <vector>
#include "polystring.h"
#include "sos.h"
#include "tests.h"

using namespace std;

int main(int argc, char * argv[]) {

    if (argc == 2 && strcmp(argv[1], "--test") == 0) {
        run_tests();
        return 0;
    }

    int d = 0;
    string positivity_condition = "PSD";
    if (argc == 2 || argc == 3) { // Input validation
        d = stoi(argv[1]);
        if (d < 0) {
            cout << "  \033[1;31m" << "Cannot use negative d value " << d << "." << "\033[0m" << endl;
            return 1;
        }
        if (argc == 3) {
            positivity_condition = argv[2];
            if (positivity_condition != "PSD" && positivity_condition != "DSOS" && positivity_condition != "SDSOS") {
                cout << "  \033[1;31m" << "Cannot use positivity condition '" << positivity_condition <<
                     "'. Can only use PSD or DSOS." << "\033[0m" << endl;
                return 1;
            }
        }
    }
    else{
        cout << "  \033[1;31m" << "Incorrect number of arguments: must be 1 or 2, not " << (argc - 1) << ".\033[0m" << endl;
        return 1;
    }

    string line, poly_string;
    vector<string> constr_strings(0, "");
    vector<string> eq_constr_strings(0, "");

    ifstream inFile;
    inFile.open("../input.txt");
    if (inFile.is_open())
    {
        int row_number = 0;
        while ( getline (inFile, line) )
            if (row_number == 0) {
                poly_string = line;  // Read objective function from first line
                row_number++;
            } else {
                constr_strings.push_back(line);  // Read constraint function into vector of constraint strings
            }
        inFile.close();
    }

//    char poly_chars[200];
//    cin.getline(poly_chars, 200, '\n');
//    string poly_string = poly_chars;
//    string poly_string = "x_1^4 - 10x_1^3 + x_1^1";
//    string poly_string = "x_n1^4 x_2^2 + x_1^2 x_2^4 - 3x_1^2 x_2^2 + 1";

    auto sol_info = sos_level_d(poly_string, constr_strings, eq_constr_strings, d, positivity_condition, 1);
    return 0;
}