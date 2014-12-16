#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include "ttmath-0.9.3\ttmath\ttmath.h"
#include <vector>

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

using namespace std;

typedef ttmath::Int<2> bigInt;

struct tpair {
	bigInt r;
	bigInt p_r;
	bigInt prime;
} ;

bigInt Mod(bigInt a, bigInt b);
bigInt Pow(bigInt a, bigInt b);
bigInt gcd(bigInt a, bigInt b);
bigInt fastModExp(bigInt b, bigInt e, bigInt m);
bigInt jacobi(bigInt a, bigInt b);
tpair tonelli(bigInt a, bigInt p);
vector<vector<bigInt>> trialDivision(bigInt n, vector<bigInt> v, int size);

#endif FUNCTIONS_H