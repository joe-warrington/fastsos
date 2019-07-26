// (c) 2019 ETH Zurich, Automatic Control Lab, Joe Warrington

#include <iostream>
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

    auto problem_tuple = read_problem_from_file("../input.txt");

    string obj_string = get<0>(problem_tuple);
    vector<string> ineq_constr_strings = get<1>(problem_tuple);
    vector<string> eq_constr_strings = get<2>(problem_tuple);
    bool read_success = get<3>(problem_tuple);

    string solver_choice = "mosek";

    if (read_success) {
        auto sol_info = sos_level_d(obj_string, ineq_constr_strings, eq_constr_strings, d, positivity_condition, 1,
                solver_choice);
        return 0;
    }
    else {
        cout << "Did not read input file successfully. Aborting.\n";
        return 1;
    }
}