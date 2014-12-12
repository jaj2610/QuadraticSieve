#include <stdio.h>
#include "functions.h"
#include "ttmath-0.9.3\ttmath\ttmath.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

using namespace std;

bigInt even = 2;

// @return a ^ b
bigInt power(bigInt a, bigInt b) 
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
	do
	{
		bigInt temp = b;
		b = mod(a, b);
		a = temp;
	} while (b != 0);

	return Abs(a);
}

// @return fast modular exponentiation of b^e mod m
bigInt exponentiate(bigInt b, bigInt e, bigInt m)
{
	bigInt n = 1;

	while (e != 0)
	{
		if (mod(e, even) == 1)
		{
			n = mod((n * b), m);
		}
		
		// e = floor(b/even)
		if (mod(b, even) == 0)
		{
			e = b / even;
		}
		else
		{
			e = (b - 1) / even;
		}

		b = mod(b * b, m);
	}

	return n;
}

// @return a Mod b
bigInt mod(bigInt a, bigInt b)
{
	return (a % b + b) % b;
}

// @return jacobi symbol (a/b)
bigInt jacobi(bigInt a, bigInt b)
{
	a = mod(a, b);
	bigInt t = 1;
	bigInt c;
	bigInt temp;
	bigInt rule1 = 8;
	bigInt rule2 = 4;

	while (a != 0) 
	{
		c = 0;
		while (mod(a, even) == 0) {
			a = a / even;
			c = bigInt(1) - c; // if c is 1, exponent is odd
		}

		if (c == 1)
		{
			if (mod(b, rule1) == 3 || mod(b, rule1) == 5)
			{
				t *= -1;
			}
		}

		if (mod(a, rule2) == 3 && mod(b, rule2) == 3) 
		{
			t *= -1;
		}

		temp = b;
		b = a;
		a = temp;
		a = mod(a, b);
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

// @return tonelli pair
tpair tonelli(bigInt a, bigInt p)
{
	bigInt b = 0;
	bigInt t = p - 1;
	bigInt s = 0;
	bigInt temp;
	tpair pair;
	bigInt c;
	bigInt i = 2;
	
	if (jacobi(a, p) == -1)
	{
		cout << "Tonelli break" << endl;
		return pair;
	}
	else
	{
		// find quad non residue
		temp = 0;
		while (temp == 1 || temp == 0)
		{
			b = b++;
			temp = jacobi(b, p);
		}

		// remove factors of 2
		while (mod(t, even) == 0)
		{
			t = t / even;
			s = s++;
		}

		c = mod(a * b * b, p);

		for (bigInt k = 1; k < s; k++)
		{
			if (exponentiate(c, power(even, (s - k - 1) * t), p) == p - 1)
			{
				i += power(even, k); 
				c = mod(c * power(b, power(even, k)), p);
			}
		}

		pair.r = mod(exponentiate(b, (i * t) / even, p) * exponentiate(a, (t + 1) / even, p), p);
		pair.p_r = p - pair.r;

		while (pair.p_r < 0)
		{
			pair.p_r + p;
		}

		pair.prime = p;
		return pair;
	}
}