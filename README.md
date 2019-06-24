# SoSchnell sum-of-squares parser

## About

SoSchnell is a parser for polynomial optimization that reads in a problem of the form 

> min f(x) = x_1 + x_2
    s.t. g_1(x) = 1 - x_1^2 - x_2^2 >= 0
    g_2(x) = -x_1x_2^2 >= 0
    
where f(x), g_1(x), ..., g_m(x) may be arbitrary polynomials in multiple variables, and solves them via a so-called *sum-of-squares relaxation*, whose degree can be specified. MOSEK (free for academic use) is the only solver supported.

## Requirements

For this code to work, the MOSEK environment variable

    MOSEKLM_LICENSE_FILE
    
on your computer must be set to the full path of your MOSEK license file. Furthermore the paths to MOSEK must also be set correctly in `CMakeLists.txt` when the code is compiled, so that the required MOSEK libraries are linked. Lastly, the environment variable

    DYLD_LIBRARY_PATH

which supplies paths to dynamically-linked libraries for the linker, must be set (or appended, if it is already set) with `<mosek install dir>/<version number>/tools/platform/osx64x86/bin` in the case of Mac OSX, and similar for other operating systems. At time of writing, only Mac OS 10.14.1 and MOSEK 9.0 have been tested together. I do not know if 

## Usage

The optimization problem is specified by entering a sequence of polynomials in the file `input.txt`. The first line is the objective function `f(x)`, and the following lines are the constraints `g_1(x)`, `g_2(x)` etc. Each constraint `g_i(x)` is interpreted as a polynomial greater than or equal to zero. This is the only kind of constraint allowed. So the file contents for the problem above would be

    x_1 + x_2
    1 - x_1^2 - x_2^2
    -x_1x_2^2
    
The variables must be called `x_1`, `x_2` etc., and only strictly positive integer exponents are allowed. Constants are also allowed. Each term of the polynomial (monomial) can be an arbitrary multiplication of `x`-variables raised to arbitrary powers, and the dimension and the degree of the polynomial will be inferred from the highest values found in the problem. Brackets and other characters are not currently interpreted.

The parser and solver are then called from the command line as follows:

    ```soschnell_cl d pos_cert```
    
In the above, `d` is an integer at least 1 that specifies the *degree* of the sum-of-squares relaxation to use. If the value of `d` is less than the minimum valid for the problem (half the highest degree of any of the polynomials, rounded up), the minimum legal value will be used, and a message to this effect displayed. 
 
The argument `pos_cert` must be either `PSD`, `SDSOS`, or `DSOS`, and specifies which positivity certificate to use throughout. These three options are in decreasing order of computational cost, but also tightness of the relaxation. They stand for *positive semidefinite*, *scaled diagonally-dominant sum-of-squares*, and *diagonally-dominant sum of squares*. 

The script then parses the polynomials described in `input.txt` and passes the resulting problem to MOSEK. The optimal value of the relaxed problem will be displayed if it was solved successfully. If the solver is not able to return an optimal solution, diagnostic information is displayed.

## Author

This code was developed by Joe Warrington at the Automatic Control Lab, ETH Zurich, Switzerland. Email <warrington@control.ee.ethz.ch>, homepage <https://people.ee.ethz.ch/~josephw>. The code is distributed under a standard GPL 3.0 license (see `LICENSE` file).