//
//  brutePaul.cpp
//  prime
//
//  Created by Paul Griffin on 18/11/14.
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
#include <thread>

using namespace std;



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

mpz_class powMod(mpz_class base, mpz_class exponent, mpz_class mod){
    mpz_class result;
    mpz_powm(result.get_mpz_t(),base.get_mpz_t(),exponent.get_mpz_t(),mod.get_mpz_t());
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
            //cout<<x<<" "<<endl;
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

auto primes = generatePrimes2(1e10);

mpz_class primalDivision(mpz_class n, vector<pair<mpz_class, long>> &factors) {
    //cout << "Primal division" << endl;
    mpz_class x = n;
    
    if (isPrime(x)) return x;

    mpz_class xroot;
    mpz_sqrt(xroot.get_mpz_t(), x.get_mpz_t());
    for (auto i: primes) {
        if (i > x) {
            break;
        }
        long count = 0;
        while (x % i == 0) {
            x /= i;
            count++;
        }
        if (count > 0) {
           // cout<<"#Primal found "<<i<<endl;
            factors.push_back(pair<mpz_class, long>(i, count));
        }
    }
    return x;
}

bool isPerfectPower(mpz_class number) {
    return mpz_perfect_power_p(number.get_mpz_t()) == 1;
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

mpz_class g(mpz_class x, mpz_class n) {
    return (x * x + 1) % n;
}

mpz_class pollard(mpz_class n, long startValue, mpz_class limit, vector<pair<mpz_class, long>> &factors) {
    mpz_class x = startValue, y = startValue, d = 1;
    auto startTime = chrono::system_clock::now();
    if (isPrime(n)) return n;
    int iterations = 0;
    while (d == 1) {
        iterations++;
        if (iterations > 10000) {
            iterations = 0;
            if (((chrono::system_clock::now() - startTime).count() / 1e6) > limit) {
                //cout << "#Hit pollard limit " << startValue << endl;
                //break;
                return n;
            }
            
        }
        x = g(x, n);
        y = g(g(y, n), n);
        mpz_class a;
        mpz_abs(a.get_mpz_t(), ((mpz_class)(x-y)).get_mpz_t());
        d = gcd(a, n);
        
        if (d != n && d != 1){
            //cout<<"#Pollard found factor "<<d<<" wtíth start value "<<startValue<<endl;
            factors.push_back(pair<mpz_class, long>(d, 1));
            n /= d;
            d = 1;
        }
    }
    return n;
}


int count2 =0;
void factorize(mpz_class base, int low,int high){
    for (int i = low; i < high; i++) {
        mpz_class number = base+i+1;
        vector<pair<mpz_class,long>> factors;
        mpz_class result = pollard(primalDivision(number, factors), 2, 10, factors);
        
        cout<<endl;
        cout<<number<<" "<<isPrime(result)<<endl;
        for ( pair<mpz_class, long> factor : factors){
            cout<<factor.first<<endl;
        }
        cout<<result<<endl;
        cout<<endl;
        if (isPrime(result)) count2++;
        
    }
}

static int c=0;
int maxNumber=10;
int parts = 8;
void threadStart(){

    int low = (maxNumber*(c+0))/parts;
    int high = (maxNumber*(c+1))/parts;
    c++;
    cout<<low<<" "<<high<<endl;
    mpz_class n("9108020935");
    mpz_class big;
    mpz_pow_ui(big.get_mpz_t(), ((mpz_class)10).get_mpz_t(), 60);
    n *= big;
    factorize(n, low, high);
}

int main(int argc, const char * argv[]) {
    

    
    thread **threads = (thread**)malloc(sizeof(thread)*8);
    
    for (int c=0; c<parts; c++){
        threads[c] = new thread(threadStart);
        chrono::milliseconds duration(100);
        this_thread::sleep_for(duration);
        
        
    }
    for (int c=0; c<parts; c++){
        threads[c]->join();
    }
    
    cout << "Factored " << count2 << endl;
    return 0;
}












