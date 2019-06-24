# SoSchnell sum-of-squares parser

## About

SoSchnell is a parser for polynomial optimization that reads in a problem of the form 

    min f(x) = x_1 + x_2
    s.t. g_1(x) = x_1^2 + x_2^2 >= 0
    g_2(x) = -x_1 + x_2^2 >= 0
    
where f(x), g_1(x), ..., g_m(x) may be arbitrary polynomials in multiple variables, and solves them via so-called *sum-of-squares relaxations*.