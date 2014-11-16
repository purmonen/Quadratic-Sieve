//
//  main.cpp
//  prime
//
//  Created by Sami Purmonen on 09/11/14.
//  Copyright (c) 2014 Sami Purmonen. All rights reserved.
//

#include <iostream>
#include <gmp.h>

#include <gmpxx.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <math.h>

using namespace std;

template <typename T>
void printVector(std::vector<T> vector) {
    for (auto t: vector) {
        std::cerr << t << " ";
    }
    std::cerr << std::endl;
}



mpz_class f(mpz_class x, mpz_class n) {
    return (x*x+1) % n;
}

mpz_class gcd(mpz_class x, mpz_class y) {
    if (x < y) {
        swap(x, y);
    }
    while (y != 0) {
        auto tmp = y;
        y = x % y;
        x = tmp;
    }
    return x;
}

bool isPrime(mpz_class x) {
    return mpz_probab_prime_p(x.get_mpz_t(), 25) == 1;
}

void writePrimesToFile(string fileName, vector<int> primes) {
    ofstream file;
    file.open(fileName);
    for (auto prime: primes) {
        file << prime << endl;
    }
}

bool isPerfectPower(mpz_class number) {
    return mpz_perfect_power_p(number.get_mpz_t()) == 1;
}

vector<int> readPrimesFromFile(string fileName) {
    vector<int> primes;
    fstream file(fileName, ios_base::in);
    int prime;
    while (file >> prime) {
        primes.push_back(prime);
    }
    return primes;
}

static vector<int> generatePrimes2(int maxprime)
{
    unsigned long long* primes;
    unsigned long long primecount=0;
    
    unsigned long long t = (maxprime>>6)+1ULL;
    primes = new unsigned long long[t];
    
    primes[0] = 0x5555555555555556ULL;
    for (int i = 1; i < ((maxprime>>6)+1ULL); i ++)
    {
        primes[i] = 0x5555555555555555ULL;
    }
    
    bool f;
    t = maxprime;
    unsigned long long maxtest = (long long)sqrtl(t) + 1;
    
    unsigned long long q;
    for ( long long i = 3; i < maxtest; i += 2)
    {
        f=primes[(i-1)>>6] & (1ULL<<(((i-1)&63)));
        if (f == true)
        {
            q = i + i;
            for( long long p = i * i-1; p < maxprime-1;p += q)
            {
                unsigned long long dc =~ (1ULL<<((p&63)));
                primes[p>>6] = primes[p>>6] & dc;
            }
        }
    }
    vector<int> result;
    
    for ( long long i = 2; i < maxprime; i ++)
    {
        f=primes[(i-1)>>6ULL] & (1ULL<<(((i-1)&63ULL)));
        if (f)
        {
            result.push_back((int)i);
            primecount++;
        }
    }
    
    return result;
}

vector<int> generatePrimes(int limit) {
    cout << "Generating primes" << endl;
    string primeFileName = "primes" + to_string(limit) + ".txt";
    
    auto primes = readPrimesFromFile(primeFileName);
    if (primes.size() > 0) {
        cout << "Done generating primes from file " << endl;
        return primes;
    }
    primes = generatePrimes2(limit);
    //    primes.push_back(2);
    //    for (int i = 3; i <= limit; i+=2) {
    //        bool isPrime = true;
    //        for (auto prime: primes) {
    //            if (i % prime == 0) {
    //                isPrime = false;
    //                break;
    //            }
    //        }
    //        if (isPrime) {
    //            primes.push_back(i);
    //        }
    //    }
    cout << "Done generating primes" << endl;
    writePrimesToFile(primeFileName, primes);
    return primes;
}

auto primes = generatePrimes(1e8);

mpz_class pollard(mpz_class n, int startValue, mpz_class limit) {
    mpz_class x = startValue, y = startValue, d = 1;
    while (d == 1) {
        if (limit-- < 0) {
            cout << "Hit pollard limit " << startValue << endl;
            return NULL;
        }
        x = f(x, n);
        y = f(f(y, n), n);
        mpz_class a;
        mpz_abs(a.get_mpz_t(), ((mpz_class)(x-y)).get_mpz_t());
        d = gcd(a, n);
    }
    return (d == n) ? NULL : d;
}

int a = 1;

mpz_class msqrtfloor(mpz_class x) {
    mpz_class xSqrt;
    mpz_class rem;
    mpz_sqrtrem(xSqrt.get_mpz_t(), rem.get_mpz_t(), x.get_mpz_t());
    return xSqrt;
}


mpz_class msqrtceiling(mpz_class x) {
    mpz_class xSqrt;
    mpz_class rem;
    mpz_sqrtrem(xSqrt.get_mpz_t(), rem.get_mpz_t(), x.get_mpz_t());
    if (rem > 0) {
        xSqrt += 1;
    }
    return xSqrt;
}


class FactorNumber {
public:
    mpz_class number;
    vector<pair<mpz_class, int>> factors;
    mpz_class quotient;
    int B;
    
    FactorNumber(mpz_class number):FactorNumber(number, factors, number) {}
    
    FactorNumber(mpz_class number, vector<pair<mpz_class, int>> factors, mpz_class quotient) {
        this->number = number;
        this->factors = factors;
        this->quotient = quotient;
        double n = number.get_d();
        B = 3*exp(0.5*sqrt(log(n)*log(log(n))));
        //        print();
    }
    
    FactorNumber pollardish(mpz_class limit) {
        cout << "Pollardish" << endl;
        vector<pair<mpz_class, int>> factors(this->factors);
        auto quotient = this->quotient;
        for (auto startValue = 2; startValue < 6; startValue++) {
            if (isPrime(quotient)) {
                cout << "IS PRIME" << endl;
                break;
            }
            auto factor = pollard(quotient, startValue, limit);
            if (factor != NULL) {
                cout << "Pollard found " << factor << endl;
                factors.push_back(std::pair<mpz_class, int>(factor, 1));
                quotient /= factor;
            }
        }
        return FactorNumber(number, factors, quotient);
    }
    
    FactorNumber trialDivision(mpz_class limit) {
        cout << "Trial division" << endl;
        mpz_class x = quotient;
        vector<pair<mpz_class, int>> factors;
        mpz_class xroot;
        mpz_sqrt(xroot.get_mpz_t(), x.get_mpz_t());
        for (auto i = 2; i <= min(xroot,limit); i++) {
            int count = 0;
            while (x % i == 0) {
                x /= i;
                count++;
            }
            if (count > 0) {
                factors.push_back(pair<mpz_class, int>(i, count));
            }
        }
        factors.insert(factors.begin(), this->factors.begin(), this->factors.end());
        return FactorNumber(number, factors, x);
    }
    
    
    FactorNumber primalDivision() {
        cout << "Primal division" << endl;
        mpz_class x = quotient;
        vector<pair<mpz_class, int>> factors;
        mpz_class xroot;
        mpz_sqrt(xroot.get_mpz_t(), x.get_mpz_t());
        for (auto i: primes) {
            if (i > x) {
                break;
            }
            int count = 0;
            while (x % i == 0) {
                x /= i;
                count++;
            }
            if (count > 0) {
                factors.push_back(pair<mpz_class, int>(i, count));
            }
        }
        factors.insert(factors.begin(), this->factors.begin(), this->factors.end());
        return FactorNumber(number, factors, x);
    }
    
    
    void print() {
        cout << "Number: " << this->number << endl;
        cout << "Factors: " << endl;
        for (auto factor: this->factors) {
            std::cout << factor.first << "^" << factor.second << "*";
        }
        cout << this->quotient << endl;
        cout << "Quotient: " << this->quotient << endl;
        cout << "B: " << this->B << endl;
        
    }
    
    void internalCheck() {
        mpz_class tmp = quotient;
        for (auto factor: factors) {
            mpz_class pow;
            mpz_pow_ui(pow.get_mpz_t(), factor.first.get_mpz_t(), factor.second);
            tmp *= pow;
        }
        if (tmp != number) {
            cout << "ERROR: is inconsistent";
        }
        
        cout << "Internal check done." << endl;
    }
    
    mpz_class Q(int x) {
        mpz_class quotientSqrt = msqrtceiling(quotient);
        quotientSqrt += x;
        quotientSqrt *= quotientSqrt;
        quotientSqrt -= quotient;
        return quotientSqrt;
    }
    
    
    mpz_class mpow(mpz_class x, int a) {
        mpz_class pow;
        mpz_pow_ui(pow.get_mpz_t(), x.get_mpz_t(), a);
        return pow;
    }
    
    FactorNumber quadraticSieve() {

        auto primeBase = generatePrimeBase();
        auto count = 20 * primeBase.size();
        vector<mpz_class> y;
        for (int x = 0; x < count; x++) {
            y.push_back(Q(x));
        }
        vector<mpz_class> oldY(y);
        vector<mpz_class> bitsets(y.size(), 0);
        for (int p = 0; p < primeBase.size(); p++) {
            auto primePair = primeBase[p];
            for (int i = primePair.second; i < y.size(); i += primePair.first) {
                y[i] /= primePair.first;
                bitsets[i] ^= ((1 << (p+1)/2));
            }
        }
        
        for (int i = 0; i < y.size(); i++) {
            if (y[i] == 1) {
                for (int p = 0; p < primeBase.size(); p+=2) {
                    cout << ((bitsets[i] >> ((p+1)/2)) & 1) << " ";
                }
                cout << endl;
                cout << "Found something good " << i << endl;
            } else {
                bitsets[i] = 0;
            }
        }
        
        for (int i = 0; i < bitsets.size(); i++) {
            for (int j = 0; j < bitsets.size(); j++) {
                for (int k = 0; k < bitsets.size(); k++) {
                    if (bitsets[i] != 0 && bitsets[j] != 0 && bitsets[k] != 0) {
                        auto zero = bitsets[i]^bitsets[j]^bitsets[k];
                        if (i != j && j != k && zero == 0) {
                            cout << "Perfect match" << endl;
                            cout << oldY[i] << " " << oldY[j] << " " << oldY[k] << endl;
                            mpz_class Y = 1;
                            for (int index: {i,j,k}) {
                                for (int p = 0; p < primeBase.size(); p += 2) {
                                    auto on = (bitsets[index] >> ((p+1)/2)) & 1;
                                    if (on == 1) {
                                        cout << primeBase[p].first << "*";
                                        Y *= primeBase[p].first;
                                    }
                                }
                            }
                            Y = msqrtceiling(Y);
                            cout << bitsets[i] << " " << bitsets[j] << " " << bitsets[k] << endl;
                            cout << "Y^2: " << Y << endl;
                            cout << i << " " << j << " " << k << endl;
                            
                            auto quotientSqrt = msqrtceiling(quotient);
                            
                            cout << i << " " << j << " " << k << endl;
                            cout << i + quotientSqrt << " " << j + quotientSqrt << " " << k + quotientSqrt << endl;
                            
                            mpz_class X = (i + quotientSqrt)*(j + quotientSqrt)*(k + quotientSqrt);
                            cout << "X^2: " << X << endl;
                            
                            cout << "N: " << quotient << endl;
                            cout << "X-Y: " << gcd((X-Y), quotient) << endl;
                            cout << "X+Y: " << gcd((X+Y), quotient) << endl;
                            
                            vector<pair<mpz_class, int>> factors(this->factors);
                            
                            mpz_class x = quotient;
                            auto commonFactors = gcd(X-Y,quotient);
                            if (commonFactors != 1) {
                                factors.push_back(pair<mpz_class, int>(commonFactors, 1));
                                x /= (commonFactors);
                            }
                            commonFactors = gcd(X+Y,quotient);
                            if (commonFactors != 1) {
                                factors.push_back(pair<mpz_class, int>(commonFactors, 1));
                                x /= (commonFactors);
                            }
                            
                            if (commonFactors != 1) {
                                return FactorNumber(number, factors, x);
                            }
                            //                            factors
                        }
                    }
                }
            }
        }
        
        return FactorNumber(number, factors, quotient);
    }
    
    vector<pair<int, int>> generatePrimeBase() {
        vector<pair<int, int>> primeBase;
        for (auto prime: primes) {
            if (prime > B) {
                break;
            }
            mpz_class tmp;
            
            mpz_mod_ui(tmp.get_mpz_t(), quotient.get_mpz_t(), prime);
            auto mod = tmp.get_si();
            for (int i = 0; i < prime; i++) {
                if ((i*i) % prime == mod) {
                    primeBase.push_back(pair<int, int>(prime, ((i - 124) % prime + prime) % prime));
                }
            }
        }
        return primeBase;
    }
};

int main(int argc, const char * argv[]) {
    mpz_class n("9011221992");
    mpz_class big;
    mpz_pow_ui(big.get_mpz_t(), ((mpz_class)10).get_mpz_t(), 60);
    n *= big;
    n += 2;
    n = 15347;

    
    vector<pair<mpz_class, int>> v;
    auto number = FactorNumber(n, v, n).quadraticSieve();
    number.internalCheck();
    number.print();
    return 0;
}