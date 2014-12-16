// Author: Jacob Jones
// Class: Cmpsc467 Fall 2014 Penn State - UP
// Date: 12/16/2014
// Description: Implementation of Quadratic Sieve for the 15 digit number: 135291536006657

#include <iostream>
#include <list>
#include <math.h>
#include "ttmath-0.9.3\ttmath\ttmath.h"
#include "functions.h"
#include <vector>

using namespace std;

extern bigInt even;
extern bigInt primes[];
extern bigInt myNumber;
extern long long myNumberLong;
extern bigInt myNumberSqrt;
extern int B;
extern int numOfPrimes;
extern int M;
extern double Threshold;
extern int factorBaseSize;
extern double thresholdDiv;

int main(int argc, const char * argv[])
{
    bigInt n = myNumber;
    bigInt sqrtN = myNumberSqrt;
    tpair p;
    bool factored = false;

    // Introduction
    cout << "Cmpsc467 Fall 2014 Penn State - UP" << endl;
    cout << "Quadratic Sieve for: " << n << "." << endl;
    cout << "Author: Jacob Jones\n" << endl;

    vector<bigInt> factorBase;
    // Calculate jacobi value for each prime less than B, and if it is 1 (n quad res mod prime) add it to base
    for (int i = 0; i < numOfPrimes; i++)
    {
        if (jacobi(n, primes[i]) == 1)
        {
            factorBase.push_back(primes[i]);
        }
    }

    // print factor base
    cout << "Factor Base: " << endl;
    for (unsigned int k = 0; k < factorBase.size(); k++)
    {
        cout << factorBase[k] << " ";
    }

    factorBaseSize = factorBase.size();

    // Tonelli's to solve congruences over the base
    vector<tpair> solved;
    for (unsigned int i = 0; i < factorBase.size(); i++)
    {
        tpair pair1 = tonelli(factorBase[i], n);
        solved.push_back(pair1);
    }

    // Print congruences
    cout << "\nTonelli Congruences: " << endl;
    for (unsigned int i = 0; i < solved.size(); i++)
    {
        cout << solved[i].prime << ": r = " << solved[i].r << ", p - r = " << solved[i].p_r << endl;
    }

    // Start the sieving process
    int sumLogSize = 2 * M;
    vector<int> sumLog;
    for (int i = 0; i < sumLogSize; i++)
    {
        sumLog.push_back(0);
    }

    bigInt temp;
    int o = 1;
    while (o < sumLogSize)
    {
        temp = sqrtN - M + o;
        for (unsigned int k = 0; k < factorBase.size(); k++)
        {
            if ((Mod(temp, factorBase[k])) == Mod(solved[k].r, factorBase[k]) || (Mod(temp, factorBase[k])) == Mod(solved[k].p_r, factorBase[k]))
            {
                sumLog[o] = sumLog[o] + (int)(.5 + log(factorBase[k].ToInt()));
            }
        }
        o++;
    }

    // Add -1 and 2 to factor base
    factorBase.insert(factorBase.begin(), 2);
    factorBase.insert(factorBase.begin(), -1);

    // reset size
    factorBaseSize = factorBase.size();

    vector<vector<bigInt>> trialDiv;
    vector<bigInt> factor1;
    vector<vector<bigInt>> factorVector;
    vector<vector<bigInt>> factorElim;

    for (int k = 1; k < sumLogSize; k++)
    {
        if (sumLog[k] >= thresholdDiv)
        {
            temp = sqrtN - M + k;
            bigInt value = Pow(temp, 2) - n;

            // Run trial division on value, and if it factors add it to factor1
            trialDiv = trialDivision(value, factorBase, factorBaseSize);
            if (!trialDiv.empty())
            {
                factor1.push_back(k);
                factorVector.push_back(trialDiv[1]);
                // Create matrix mod 2
                vector<bigInt> modVect;
                for (int l = 0; l < factorBaseSize; l++)
                {
                    modVect.push_back(Mod(trialDiv[1].at(l), even));
                }
                factorElim.push_back(modVect);
            }
        }
    }

    // Gaussian Elim
    // add identity matrix
    for (unsigned int l = 0; l < factor1.size(); l++)
    {
        for (unsigned int k = 0; k < factor1.size(); k++)
        {
            if (l == k)
            {
                factorElim[l].push_back(1);
            }
            else
            {
                factorElim[l].push_back(0);
            }
        }
    }

    // Perform Elim
    for (int col = 0; col < factorBaseSize; col++)
    {
        for (unsigned int row = 0; row < factorElim.size(); row++)
        {
            if (factorElim[row].at(col) == 1)
            {
                vector<bigInt> sum = factorElim[row];
                for (unsigned int row2 = row + 1; row2 < factorElim.size(); row2++)
                {
                    if (factorElim[row2].at(col) == 1)
                    {
                        vector<bigInt> sum2 = factorElim[row2];
                        vector<bigInt> output;
                        for (unsigned int w = 0; w < sum.size(); w++)
                        {
                            output.push_back(Mod((sum2[w] + sum[w]), even));
                        }
                        factorElim[row2].clear();
                        for (unsigned int col2 = 0; col2 < output.size(); col2++)
                        {
                            factorElim[row2].push_back(output[col2]);
                        }
                    }
                    // Remove initial row
                    factorElim[row].clear();
                    factorElim.erase(factorElim.begin() + row);
                    break;
                }
            }
        }
    }

    // Solve x^2 congruent to y^2 (mod n)
    for (unsigned int i = 0; i < factorElim.size(); i++)
    {
        bigInt lhs = 1;
        bigInt rhs = 1;
        bigInt x;
        bigInt y;
        vector<bigInt> sum;

        // initialize sum to 0
        for (int k = 0; k < factorBaseSize; k++)
        {
            sum.push_back(0);
        }

        // Calculate x value by taking product of rows
        for (unsigned int j = factorBaseSize; j < factorElim[0].size(); j++)
        {
            if (factorElim[i].at(j) == 1)
            {
                int pos = j - factorBaseSize;
                bigInt value = factor1[pos];
                bigInt xValue = sqrtN - M + value;
                lhs = lhs * xValue;
                lhs = Mod(lhs, n);

                vector<bigInt> tempVect1 = factorVector[pos];
                vector<bigInt> tempVect2;
                for (int w = 0; w < factorBaseSize; w++)
                {
                    tempVect2.push_back(sum[w] + tempVect1[w]);
                }
                for (int u = 0; u < factorBaseSize; u++)
                {
                    sum[u] = tempVect2[u];
                }
            }
            vector<bigInt> finalSum;
            for (int m = 0; m < factorBaseSize; m++)
            {
                finalSum.push_back(sum[m] / even);
            }

            for (int m = 0; m < factorBaseSize; m++)
            {
                y = fastModExp(factorBase[m], finalSum[m], n);
                rhs = rhs * y;
                rhs = Mod(rhs, n);
            }

            long long gcdValue = (lhs - rhs).ToInt();

            // find gcd(x-y, n)
            long long GCD = gcd(gcdValue, myNumberLong);

            if ((GCD > 1) && (GCD < myNumberLong))
            {
                factored = true;
                cout << "Factored successfully." << endl;
                cout << GCD << " * " << n / GCD << endl;
                cout << "M: " << M << ", B: " << B << ", factor base size: " << factorBaseSize << endl;
                break;
            }
        }

        if (!factored)
        {
            cout << "Failed to factor..." << endl;
        }
    }


     system("pause");

    return 0;
}