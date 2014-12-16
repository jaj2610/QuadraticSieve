#include <stdio.h>
#include "functions.h"
#include "ttmath-0.9.3\ttmath\ttmath.h"
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

bigInt even = 2;
bigInt primes[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149,
151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,
317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499,
503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691,
701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907,
911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997 };

bigInt myNumber = 135291536006657; // given by assignment
bigInt myNumberSqrt = 11631488; // floor(sqrt(myNumber))... sqrt(myNumber) = 11631488.9849...
int B = 205; // calculated via formula
int numOfPrimes = 40; // pi(B)
int M = 6000; // calculated based on B
double Threshold = 1.5; // copied from book methodology
int factorBaseSize;
double thresholdDiv = (.5 * log(135291536006657)) + log(M) - (Threshold * log(B));

// @return a ^ b
bigInt Pow(bigInt a, bigInt b) 
{
	bigInt temp = 1;
	for (bigInt k = 0; k < b; k++) 
	{
		temp *= a;
	}
	return temp;
}

// @return gcd(a, b)
bigInt gcd(bigInt a, bigInt b)
{
	while (b != 0)
	{
		bigInt temp = b;
		b = Mod(a, b);
		a = temp;
	}

	return Abs(a);
}

// @return fast modular exponentiation of b^e mod m
bigInt fastModExp(bigInt b, bigInt e, bigInt m)
{
	bigInt n = 1;

	while (e != 0)
	{
		if (Mod(e, even) == 1)
		{
			n = Mod((n * b), m);
		}
		
		// e = floor(b/even)
		if (Mod(b, even) == 0)
		{
			e = b / even;
		}
		else
		{
			e = (b - 1) / even;
		}

		b = Mod(b * b, m);
	}

	return n;
}

// ensure modulus doesn't return negative
// @return a Mod b
bigInt Mod(bigInt a, bigInt b)
{
	bigInt output = a % b;
	if (output < 0)
	{
		output = output + b;
	}

	return output;
}

// @return jacobi symbol (a/b)
bigInt jacobi(bigInt a, bigInt b)
{
	a = Mod(a, b);
	bigInt t = 1;

	while (a != 0)
	{
		bigInt c = 0;
		while (Mod(a, even) == 0)
		{
			a = a / 2;
			c = (bigInt)1 - c;
		}
		if (c == 1)
		{
			if (Mod(b, (bigInt)8) == 3 || (Mod(b, (bigInt)8) == 5))
			{
				t = (bigInt)-1 * t;
			}
		}
		if ((Mod(a, (bigInt)4) == 3) && (Mod(b, (bigInt)4) == 3))
		{
			t = (bigInt)-1 * t;
		}
		bigInt temp = b;
		b = a;
		a = temp;
		a = Mod(a, b);
	}
	if (b == 1)
	{
		return t;
	}
	else
	{
		return 0;
	}
}

// @param p odd prime
// @param a integer 1 < a < p
// @return tonelli pair
tpair tonelli(bigInt p, bigInt a)
{
	if (a > (p - 1))
	{
		a = Mod(a, p);
	}

	if (((jacobi(a, p)) == -1) || ((jacobi(a, p)) == 0))
	{
		cout << a << " is a quad residue mod " << p << endl;
	}

	bigInt b;
	for (b = 0; b < p; b++)
	{
		if ((jacobi(b, p)) == -1)
		{
			break;
		}
	}

	bigInt s = 0;
	bigInt t = p - 1;
	while (Mod(t, even) == 0)
	{
		t = t / even;
		s = s + 1;
	}

	bigInt i = 2;
	bigInt c = Mod(Mod(a, p) * Mod(b, p) * Mod(b, p), p);

	for (bigInt k = 1; k < s; k++)
	{
		bigInt e = (bigInt)(t * Pow(even, s - k - 1));
		if (fastModExp(c, e, p) == -1 || fastModExp(c, e, p) == (p - 1))
		{
			i = (i + (bigInt)Pow(even, k));
			c = Mod(Mod(c, p) * fastModExp(b, (bigInt)Pow(even, k), p), p);
		}
	}

	tpair pair;
	bigInt temp1 = fastModExp(b, (bigInt)((i * t) / even), p);
	bigInt temp2 = (bigInt) fastModExp((bigInt)a, (bigInt)((t + 1) / even), (bigInt)p);
	pair.r = (bigInt)Mod((temp1 * temp2), p);
	pair.p_r = p - pair.r;

	return pair;
}

// @return exponent vector
vector<vector<bigInt>> trialDivision(bigInt n, vector<bigInt> v, int size)
{
	bigInt f = n;
	vector<bigInt> p;
	vector<bigInt> e;

	// Init p to be factor base
	for (int i = 0; i < size; i++)
	{
		p[i] = v[i];
		e[i] = 0;
	}

	if (f < 0)
	{
		e[0] = 1;
		f = Abs(f);
	}

	for (int i = 1; i < size; i++)
	{
		bigInt d = p[i];
		while (Mod(f, d) == 0)
		{
			f = f / d;
			e[i] = e[i] + 1;
		}
	}

	vector<vector<bigInt>> output;

	output.push_back(p);
	output.push_back(e);
	
	return output;
}