# FastSOS sum-of-squares parser

## About

FastSOS is a parser for polynomial optimization that reads in a problem of the form 

> min f(x) = x1 + x2,

> s.t. g_1(x) = 1 - x1^2 - x2^2 >= 0, g_2(x) = -x1x2^2 >= 0

> h_1(x) = x1 - x1^2 = 0
    
where f(x), g_1(x), ..., g_m(x), h_1(x), ..., h_p(x) may be arbitrary polynomials in multiple variables, and solves them via a so-called *sum-of-squares relaxation*, whose degree can be specified. MOSEK (free for academic use) is the only solver supported.

## Requirements

### Compiling and linking

When compiling the code, the paths to MOSEK must also be set correctly in `CMakeLists.txt` when the code is compiled, so that the required MOSEK libraries are linked. The environment variable

    DYLD_LIBRARY_PATH

which supplies paths to dynamically-linked libraries for the linker, must also be set (or appended, if the variable is already in use) with `<mosek install dir>/<version number>/tools/platform/osx64x86/bin` in the case of Mac OSX, and similar for other operating systems. At time of writing, the only compilation tested has been on Mac OS 10.14.1 with MOSEK 9.0.

### Running

For this code to run, the MOSEK environment variable

    MOSEKLM_LICENSE_FILE

must be set to the full path of your MOSEK license file. On Mac OSX, the two environment variables mentioned above can be set up for permanent use in all command line sessions by creating or editing `~/.bash_profile`.

## Usage

### Input format

The optimization problem is specified by entering a sequence of polynomials in the file `input.txt`. The objective function, inequality and equality constraints are indicated by an initial `f ` (f plus a space), `g ` or `h ` respectively. There must be exactly one objective line. The inequality constraints g_i(x) are interpreted as _greater than or equal to_ zero, in accordance with much of the polynomial optimization literature. Lines not starting with the correct letter-plus-space combo are ignored. So the example problem above would be entered as follows:

    f x1 + x2
    g 1 - x1^2 - x2^2
    g -x1x2^2
    h x1 - x1^2
    
The variables must be called `x1`, `x2` etc., and only strictly positive integer exponents are allowed. Constants are also allowed. Each term of the polynomial (monomial) can be an arbitrary multiplication of `x`-variables raised to arbitrary powers, and the dimension and the degree of the polynomial will be inferred from the highest values found in the problem. Brackets and other characters are not currently interpreted. Unused coordinates, for example `x2` not appearing in a problem where `x1` and `x3` appear, are automatically eliminated to reduce computational load.

### Optimization

The parser and solver are then called from the command line as follows:

```
fastsos d pos_cert
```
    
In the above, `d` is an integer at least 1 that specifies the *degree* of the sum-of-squares relaxation to use. If the value of `d` is less than the minimum valid for the problem (half the highest degree amongst all the polynomials in the problem, rounded up), the minimum legal value will be used, and a message to this effect displayed. 
 
The argument `pos_cert` must be either `PSD`, `SDSOS`, or `DSOS`, and specifies which positivity certificate to use throughout. These three options are in decreasing order of computational cost, but also generate weaker and weaker relaxations. They stand for *positive semidefinite*, *scaled diagonally-dominant sum-of-squares*, and *diagonally-dominant sum of squares*. Strictly the latter two should simply be called (scaled) diagonally dominant or (S)DD, as they refer to matrix positivity conditions. However (S)DSOS has become the more common terminology.

The script then parses the polynomials described in `input.txt` and passes the resulting problem to MOSEK. The optimal value of the relaxed problem will be displayed if it was solved successfully. If the solver is not able to return an optimal solution, diagnostic information is displayed.

### Tests

The polynomial parser and simple SOS programming behaviour can be tested using the following command: 

```
fastsos --test
``` 

This checks the results of parsing certain and solving certain hard-coded SOS problems against expected behaviour, and reports on any failures.

### Interfaces

The intention at time of writing is to incorporate the code into a Python library for ease of access in other high-level scripts. However this has not yet been implemented.

## Author

This code was developed by Joe Warrington at the Automatic Control Lab, ETH Zurich, Switzerland. Email <warrington@control.ee.ethz.ch>, homepage <https://people.ee.ethz.ch/~josephw>. The code is distributed under a standard GNU General Public License v3.0; see `LICENSE` file.