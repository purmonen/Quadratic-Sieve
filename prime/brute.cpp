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
#include <assert.h>
#include <time.h>

using namespace std;

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

mpz_class powMod(mpz_class base, mpz_class exponent, mpz_class mod){
    mpz_class result;
    mpz_powm(result.get_mpz_t(),base.get_mpz_t(),exponent.get_mpz_t(),mod.get_mpz_t());
    return result;
}

mpz_class powNoMod(mpz_class base, long exponent){
    mpz_class result;
    mpz_pow_ui(result.get_mpz_t(),base.get_mpz_t(),exponent);
    return result;
}

static gmp_randclass gmpRandom(gmp_randinit_default);
bool millerRabin(mpz_class n, int tries){
    if (n <= 1) return false;
    if (n == 2) return true;
    if (n == 3) return true;
    if (n%2 == 0) return false;
    mpz_class s = 0;
    mpz_class d = n-1;
    while ( d%2 == 0){
        d /= 2;
        s++;
    }
    for (int i=0; i<tries; i++){
        mpz_class a = gmpRandom.get_z_range(n-2)+2;
        mpz_class x = powMod(a, d, n);
        if (x == 1 || x == n-1) continue;
        bool ok = false;
        for (int j = 1; j<s; j++){
            x = (x * x) % n;
            cout<<x<<" "<<endl;
            if (x == 1) return false;
            if (x == (n-1)) {
                ok = true;
                break;
            }
        }
        if (!ok) return false;
    }
    
    
    return true;
}

bool isPrime(mpz_class x) {
    bool a = mpz_probab_prime_p(x.get_mpz_t(), 25) >= 1;
    bool b = millerRabin(x, 25);
    if (a != b) {
        cout << "x " << x << endl;
    }
    assert(a==b);
    return a;
    
}

bool isPerfectPower(mpz_class number) {
    return mpz_perfect_power_p(number.get_mpz_t()) == 1;
}

static vector<long> generatePrimes2(long maxprime)
{
    unsigned long long* primes;
    unsigned long long primecount=0;
    
    unsigned long long t = (maxprime>>6)+1ULL;
    primes = new unsigned long long[t];
    
    primes[0] = 0x5555555555555556ULL;
    for (long i = 1; i < ((maxprime>>6)+1ULL); i ++)
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
    vector<long> result;
    
    for ( long long i = 2; i < maxprime; i ++)
    {
        f=primes[(i-1)>>6ULL] & (1ULL<<(((i-1)&63ULL)));
        if (f)
        {
            result.push_back((long)i);
            primecount++;
        }
    }
    
    return result;
}


vector<long> generatePrimes(long limit) {
    auto primes = generatePrimes2(limit);
    return primes;
}

auto primes = generatePrimes(1e8);

class PollardInfo {
public:
    mpz_class x;
    mpz_class y;
    mpz_class n;
    mpz_class limit;
    mpz_class factor;
    bool overLimit;
    
    PollardInfo(mpz_class x, mpz_class y, mpz_class n, mpz_class limit, mpz_class factor, bool overLimit) {
        this->x = x,
        this->x = y;
        this->n = n;
        this->limit = limit;
        this->factor = factor;
        this->overLimit = overLimit;
    }
};

mpz_class f(mpz_class x, mpz_class n) {
    return (x*x+1) % n;
}

vector<pair<mpz_class, long>> primalDivision(mpz_class &quotient) {
    vector<pair<mpz_class, long>> factors;
    mpz_class x = quotient;
    mpz_class xroot;
    mpz_sqrt(xroot.get_mpz_t(), x.get_mpz_t());
    for (auto i: primes) {
        if (i > xroot) {
            break;
        }
        long count = 0;
        while (x % i == 0) {
            x /= i;
            count++;
        }
        if (count > 0) {
            factors.push_back(pair<mpz_class, long>(i, count));
        }
    }
    factors.insert(factors.begin(), factors.begin(), factors.end());
    
    quotient = x;
    return factors;
}


PollardInfo pollard(PollardInfo pollardInfo) {
    mpz_class x = pollardInfo.x, y = pollardInfo.y, n = pollardInfo.n; mpz_class d = 1;
    auto startTime = chrono::system_clock::now();
    int iterations = 0;
    while (d == 1) {
        iterations++;
        if (iterations > 10000) {
            iterations = 0;
//            cout << "Time " << (chrono::system_clock::now() - startTime).count() << endl;
            if (((chrono::system_clock::now() - startTime).count() / 1e6) > pollardInfo.limit) {
                cout << "Hit pollard limit " << endl;
                pollardInfo.x = x;
                pollardInfo.y = y;
                pollardInfo.overLimit = true;
                return pollardInfo;
            }
        }
        x = f(x, n);
        y = f(f(y, n), n);
        mpz_class a;
        mpz_abs(a.get_mpz_t(), ((mpz_class)(x-y)).get_mpz_t());
        d = gcd(a, n);
    }
    
    pollardInfo.x = x;
    pollardInfo.y = y;
    pollardInfo.factor = (d == n) ? NULL : d;
    return pollardInfo;
}

class Logger {
public:
    std::chrono::time_point<std::chrono::system_clock> lastTime = chrono::system_clock::now();
    bool isFirst = true;
    string lastLog = "";
    void log(string logString) {
        auto currentTime = chrono::system_clock::now();
        auto duration = (currentTime - lastTime) / 1e6;
        lastTime = currentTime;
        if (!isFirst) {
            cout << lastLog << "\t" << duration.count() << "s" << endl;
        }
        isFirst = false;
        lastLog = logString;
        cout << logString << endl;
    }
};

int main(int argc, const char * argv[]) {
    cout << "Starting factoring" << endl;
    mpz_class startN("9011221992");
    mpz_class big;
    mpz_pow_ui(big.get_mpz_t(), ((mpz_class)10).get_mpz_t(), 60);
    startN *= big;
    //    n = 15348;
    int count = 0;
    
    mpz_class startValue = 3;
    
    int totalNumbers = 10;
    for (int i = 1; i < totalNumbers; i++) {
        mpz_class n = startN + i;
        primalDivision(n);
        auto pollardInfo = PollardInfo(startValue, startValue, n, 1000, NULL, false);
        if (n == 1) {
            cout << "Found by primal division" << endl;
            count++;
            continue;
        }
        cout << "Factoring " << n << endl;
        
        pollardInfo = pollard(pollardInfo);
        
        if (pollardInfo.factor != NULL) {
            count++;
            cout << "Found factor " << pollardInfo.factor << endl;
            if (isPrime(pollardInfo.factor)) {
                cout << "IS PRIME" << endl;
            }
            pollardInfo.factor = NULL;
        } else {
//            cout << "Pollard failed " << pollardInfo.factor << endl;
        }
    }
    cout << "Factored " << count << "/" << totalNumbers << endl;
    return 0;
}