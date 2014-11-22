//
//  test.cpp
//  prime
//
//  Created by Paul Griffin on 20/11/14.
//  Copyright (c) 2014 Sami Purmonen. All rights reserved.
//

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
#include <fstream>
#include <thread>
#include <chrono>
#include <algorithm>

using namespace std;

template <class T>
T gcd(T x, T y) {
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

static long modPow(long base,long power, long mod){
    long result = base;
    while (power>1) {
        result *= base;
        result %= mod;
    }
    return result;
}


int main(int argc, const char * argv[]) {
    
    long max = 1e6;
    
    auto primes = generatePrimes2(max);
    
    
    for (long n = 3; n < max; n++){
        bool isCarmichael = true;
        for (auto p:primes){
            if (p > n) break;
            if (((n-1)%(p-1) != 0 && n%p == 0)|| p==n ){
                isCarmichael = false;
                break;
            }
        }
        if (isCarmichael){
            long value = n;
            for (auto p:primes)
                if (value%p==0 && value/p%p ==0) {isCarmichael = false; break;};
            if (isCarmichael) cout<<n<<" ";
        }
    }
    cout<<endl;
    return 1;
}














