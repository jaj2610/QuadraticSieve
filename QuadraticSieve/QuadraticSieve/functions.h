#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include "ttmath-0.9.3\ttmath\ttmath.h"

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

typedef ttmath::Int<2> bigInt;

struct tpair {
	bigInt r;
	bigInt p_r;
	bigInt prime;
};

bigInt mod(bigInt a, bigInt b);
bigInt power(bigInt a, bigInt b);
bigInt gcd(bigInt a, bigInt b);
bigInt exponentiate(bigInt b, bigInt e, bigInt m);
bigInt jacobi(bigInt a, bigInt b);
tpair tonelli(bigInt a, bigInt p);

#endif FUNCTIONS_H