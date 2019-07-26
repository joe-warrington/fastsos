//
// Created by Joe Warrington on 2019-07-26.
//

#include <iostream>
#include <iomanip>
#include <sstream>
#include <sys/time.h>
#include <math.h>
#include "timing.h"

using namespace std;

typedef unsigned long long timestamp_t;
timestamp_t timenow() {
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

string time_string(unsigned long long us_in) {
    ostringstream time_string_stream;     // Create an output string stream
    time_string_stream << std::fixed;    // Set Fixed -Point Notation
    time_string_stream << setprecision(3);    // Set precision to 2 digits

//    if (log10((double) us_in) < 3) {
//        time_string_stream << us_in;
//        return time_string_stream.str() + " us";
//    }
//    else
    if (log10((double) us_in) >= 6) {
        time_string_stream << us_in / 1000000.0;
        return time_string_stream.str() + " s";
    }
    else {
        time_string_stream << us_in / 1000.0;
        return time_string_stream.str() + " ms";
    }
}