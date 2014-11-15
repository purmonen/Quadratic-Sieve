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

using namespace std;

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

vector<int> generatePrimes(int limit) {
    cout << "Generating primes" << endl;
    vector<int> primes;
    primes.push_back(2);
    for (int i = 3; i <= limit; i+=2) {
        bool isPrime = true;
        for (auto prime: primes) {
            if (i % prime == 0) {
                isPrime = false;
                break;
            }
        }
        if (isPrime) {
            primes.push_back(i);
        }
    }
    cout << "Done generating primes" << endl;
    
    return primes;
}

auto primes = generatePrimes(1e5);

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

class FactorNumber {
public:
    mpz_class number;
    vector<pair<mpz_class, int>> factors;
    mpz_class quotient;
    
    
    FactorNumber(mpz_class number) {
        this->number = number;
        this->quotient = number;
    }
    
    FactorNumber(mpz_class number, vector<pair<mpz_class, int>> factors, mpz_class quotient) {
        this->number = number;
        this->factors = factors;
        this->quotient = quotient;
        print();
    }
    
    FactorNumber pollardish(mpz_class limit) {
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
        mpz_class x = quotient;
        vector<pair<mpz_class, int>> factors;
        mpz_class xroot;
        mpz_sqrt(xroot.get_mpz_t(), x.get_mpz_t());
        for (auto i: primes) {
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
        
    }
};

int main(int argc, const char * argv[]) {
    mpz_class n("9011221992");
    mpz_class big;
    mpz_pow_ui(big.get_mpz_t(), ((mpz_class)10).get_mpz_t(), 60);
    n *= big;
    n += 1;
    
    auto number = FactorNumber(n).primalDivision().pollardish(1e5);
    return 0;
}