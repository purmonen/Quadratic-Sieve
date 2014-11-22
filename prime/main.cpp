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


#define GMP_LIMB_BITS_SHIFTER 6
struct bitarray{
    mp_limb_t* limb;
    long bitSize;
    mp_size_t limbCount;
    
    bitarray(long bitCount){
        bitSize = bitCount;
        //assert(bitCount!=0);
        //cout<<"creating with normal constructor"<<endl;
        limbCount = (bitCount>>GMP_LIMB_BITS_SHIFTER)+1;
        //cout<<"allocating "<<limbCount<<" for "<<bitCount<<"bits"<<endl;
        limb = new mp_limb_t[limbCount];
        memset(limb, 0, limbCount*GMP_LIMB_BITS>>3);
    }
    
    bitarray(bitarray const &b){
       // cout<<"creating with copy constructor"<<endl;
        limbCount = b.limbCount;
        bitSize = b.bitSize;
        limb = new mp_limb_t[limbCount];
        memcpy(limb, b.limb, limbCount*8);
    }
    
    ~bitarray(){
        //cout<<"destorying"<<endl;
        delete[] limb;
    }
    
    
};
void swap(bitarray &a, bitarray &b){
    swap(a.limb,b.limb);
    swap(a.limbCount,b.limbCount);
}

void setBit(bitarray &x, const long &i) {
    //assert(i<x.bitSize);
    const long limb = i >> GMP_LIMB_BITS_SHIFTER;
    const long bitInLimb = i & 0b11111111ull;
    x.limb[limb] |= 1ull << bitInLimb;
}

void unsetBit(bitarray &x, const long &i) {
    //assert(i<x.bitSize);
    const long limb = i >> GMP_LIMB_BITS_SHIFTER;
    const long bitInLimb = i & 0b11111111ull;
    x.limb[limb] &= ~(1ull << bitInLimb);
}

void flipBit(bitarray &x, const long &i) {
    //assert(i<x.bitSize);
    const long limb = i >> GMP_LIMB_BITS_SHIFTER;
    const long bitInLimb = i & 0b11111111ull;
    x.limb[limb] ^= (1ull << bitInLimb);
}

inline const bool isBitSet(const bitarray &x, const long &i) {
    //assert(i<x.bitSize);
    const long limb = i >> GMP_LIMB_BITS_SHIFTER;
    const long bitInLimb = i & 0b11111111ull;
    return (x.limb[limb]>>bitInLimb) & 1;
}

long countOnes(bitarray &x) {
    return mpn_popcount(x.limb,x.limbCount);
}

bool isZero(bitarray &bitset) {
    auto i = bitset.limbCount-1;
    for (; i>=0; i--) {
        if (bitset.limb[i] != 0) return true;
    }
    return false;}

long leftMostOne(bitarray &bitset) {
    auto i = bitset.limbCount-1;
    for (; i>=0; i--) {
        if (bitset.limb[i] != 0) break;
    }
    //cout<<"found in limb "<<i<<endl;
    auto x = bitset.limb[i];
    long count = -1;
    while (x != 0) {
        count++;
        x = x >> 1;
    }
    return (i<<GMP_LIMB_BITS_SHIFTER) + count;
}

template <typename T>
void printVector(std::vector<T> vector) {
    for (auto t: vector) {
        std::cerr << t << " ";
    }
    std::cerr << std::endl;
}

void printBitVector(std::vector<bitarray> vector, size_t length) {
    
    for (auto v: vector) {
        for (int i = 0; i < length; i++) {
            cout << isBitSet(v, i);
        }
        cout << endl;
    }
    cout << std::endl;
}

void printBits(bitarray v, int length =-1) {
    if (length==-1)
        length = v.limbCount*64;
    for (int i = 0; i < length; i++) {
        cout << isBitSet(v, i);
    }
    cout << std::endl;
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

mpz_class powMod(const mpz_class &base,const mpz_class &exponent,const mpz_class &mod){
    mpz_class result;
    mpz_powm(result.get_mpz_t(),base.get_mpz_t(),exponent.get_mpz_t(),mod.get_mpz_t());
    return result;
}

mpz_class powNoMod(const mpz_class &base, const long &exponent){
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



void writePrimesToFile(string fileName, vector<long> primes) {
    ofstream file;
    file.open(fileName);
    for (auto prime: primes) {
        file << prime << "\n";
    }
}

bool isPerfectPower(mpz_class number) {
    return mpz_perfect_power_p(number.get_mpz_t()) == 1;
}

mpz_class legendreSymbol(const mpz_class &n, const mpz_class &p){
    mpz_class tmp = (p-1)/2;
    mpz_powm(tmp.get_mpz_t(), n.get_mpz_t(), tmp.get_mpz_t(), p.get_mpz_t());
    return tmp > 1 ? -1 : tmp;
}

vector<long> readPrimesFromFile(string fileName) {
    vector<long> primes;
    fstream file(fileName, ios_base::in);
    long prime;
    while (file >> prime) {
        primes.push_back(prime);
    }
    return primes;
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
    cout << "Generating primes" << endl;
    string primeFileName = "primes" + to_string(limit) + ".txt";
    
    auto primes = readPrimesFromFile(primeFileName);
    if (primes.size() > 0) {
        cout << "Done generating primes from file " << endl;
        return primes;
    }
    primes = generatePrimes2(limit);
    //    primes.push_back(2);
    //    for (long i = 3; i <= limit; i+=2) {
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
    //    writePrimesToFile(primeFileName, primes);
    return primes;
}

auto primes = generatePrimes(1e7);

mpz_class pollard2(mpz_class n, long startValue, mpz_class limit, vector<pair<mpz_class, long>> &factors) {
    mpz_class x = startValue, y = startValue, d = 1;
    auto startTime = chrono::system_clock::now();
    if (isPrime(n )) return n;
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
        x = f(x, n);
        y = f(f(y, n), n);
        mpz_class a;
        mpz_abs(a.get_mpz_t(), ((mpz_class)(x-y)).get_mpz_t());
        d = gcd(a, n);
        
        if (d != n && d != 1){
            //cout<<"#Pollard found factor "<<d<<" wtÃ­th start value "<<startValue<<endl;
            factors.push_back(pair<mpz_class, long>(d, 1));
            n /= d;
            if (isPrime( n)) return n;
            d = 1;
        }
    }
    return n;
}

mpz_class pollard(const mpz_class &n, const long &startValue,const  mpz_class &limit) {
    mpz_class x = startValue, y = startValue, d = 1;
    auto startTime = chrono::system_clock::now();
    int iterations = 0;
    while (d == 1) {
        iterations++;
        if (iterations > 1e7) {
            iterations = 0;
            if (((chrono::system_clock::now() - startTime).count() / 1e6) > limit) {
                cout << "Hit pollard limit " << startValue << endl;
                return NULL;
            }
        }
        x = f(x, n);
        y = f(f(y, n), n);
        mpz_class a;
        mpz_abs(a.get_mpz_t(), ((mpz_class)(x-y)).get_mpz_t());
        d = gcd(a, n);
    }
    return (d == n) ? NULL : d;
}



//long a = 1;

mpz_class msqrtfloor(const mpz_class &x) {
    mpz_class xSqrt;
    mpz_class rem;
    mpz_sqrtrem(xSqrt.get_mpz_t(), rem.get_mpz_t(), x.get_mpz_t());
    return xSqrt;
}


mpz_class msqrtceiling(const mpz_class &x) {
    mpz_class xSqrt;
    mpz_class rem;
    mpz_sqrtrem(xSqrt.get_mpz_t(), rem.get_mpz_t(), x.get_mpz_t());
    if (rem > 0) {
        xSqrt += 1;
    }
    return xSqrt;
}

bool hasRoot(mpz_class x, long p){
    mpz_class tmp;
    mpz_powm_ui(tmp.get_mpz_t(),x.get_mpz_t(),(p-1)/2,mpz_class(p).get_mpz_t());
    //        mpz_pow_ui(tmp.get_mpz_t(),x.get_mpz_t(),(p-1)/2);
    //        mpz_mod_ui(tmp.get_mpz_t(),tmp.get_mpz_t(),p);
    return tmp == 1;
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

class FactorNumber {
public:
    mpz_class number;
    vector<pair<mpz_class, long>> factors;
    mpz_class quotient;
    mpz_class quotientSqrt;
    long B;
    
    FactorNumber(mpz_class number):FactorNumber(number, factors, number) {}
    
    FactorNumber(mpz_class number, vector<pair<mpz_class, long>> factors, mpz_class quotient) {
        this->number = number;
        this->factors = factors;
        this->quotient = quotient;
        this->quotientSqrt = msqrtceiling(quotient);
        double n = quotient.get_d();
        B = 1.8 * exp(0.5*sqrt(log(n)*log(log(n)))) + 500;
        //        print();
        
        if (isPrime(quotient)) {
            cout << "Quotient is prime" << endl;
            this->factors.push_back(pair<mpz_class, long>(quotient, 1));
            this->quotient = 1;
        }
    }
    
    
    FactorNumber pollardish(mpz_class limit,int startvalue = 2) {
        
        cout << "Pollardish" << endl;
        vector<pair<mpz_class, long>> factors(this->factors);
        auto quotient = this->quotient;
        
        
        if (!isPrime(this->quotient) && this->quotient != 1) {
            for (auto startValue = startvalue; startValue <=4 ; startValue++) {
                if (isPrime(quotient)) {
                    cout << "IS PRIME" << endl;
                    break;
                }
                auto factor = pollard(quotient, startValue, limit);
                if (factor != NULL) {
                    cout << "Pollard found " << factor << endl;
                    factors.push_back(std::pair<mpz_class, long>(factor, 1));
                    quotient /= factor;
                    return FactorNumber(number, factors, quotient).pollardish(limit,startValue);
                }
            }
        }
        return FactorNumber(number, factors, quotient);
    }
    
    FactorNumber trialDivision(mpz_class limit) {
        cout << "Trial division" << endl;
        mpz_class x = quotient;
        vector<pair<mpz_class, long>> factors;
        mpz_class xroot;
        mpz_sqrt(xroot.get_mpz_t(), x.get_mpz_t());
        for (auto i = 2; i <= min(xroot,limit); i++) {
            long count = 0;
            while (x % i == 0) {
                x /= i;
                count++;
            }
            if (count > 0) {
                factors.push_back(pair<mpz_class, long>(i, count));
            }
        }
        factors.insert(factors.begin(), this->factors.begin(), this->factors.end());
        return FactorNumber(number, factors, x);
    }
    
    
    FactorNumber primalDivision(long maxValue = 1e13) {
        cout << "Primal division" << endl;
        mpz_class x = quotient;
        vector<pair<mpz_class, long>> factors;
        mpz_class xroot;
        mpz_sqrt(xroot.get_mpz_t(), x.get_mpz_t());
        for (auto i: primes) {
            if (i>maxValue) break;
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
    
    mpz_class mpow(const mpz_class &x, const long &a) {
        mpz_class pow;
        mpz_pow_ui(pow.get_mpz_t(), x.get_mpz_t(), a);
        return pow;
    }
    
    inline const mpz_class multiPolynomial(const mpz_class &x) {
        return (quotientSqrt + x)*(quotientSqrt + x) - quotient;
    }
    
    
    void sieve(
               //Immutable parameters
               const long lowLimit,
               const long highLimit,
               const vector<float> &oldPrimeLogs,
               const vector<pair<long, vector<long>>> &primeBase,
               const vector<long> &oldPrime,
               const vector<long> &oldRoot,
               const long maxIndex,
               const long bitWidth,
               
               //Mutable parameters
               long oldYLogs[],
               vector<mpz_class> &oldY,
               vector<mpz_class> &oldX,
               long &sieveCount,
               vector<bitarray> &bitsets,
               std::atomic_flag &lock,
               vector<long> &oldi){
        
       
        
        
        

        
        const auto maxPrimeLog = oldPrimeLogs[oldPrimeLogs.size()-1];
        auto nextIndex = lowLimit;
        long lastLog = 0;
        for (auto x= lowLimit;x<highLimit;x++){
            if (x >= nextIndex) {
                mpz_class value = multiPolynomial(x);
                lastLog = mpz_sizeinbase(value.get_mpz_t(), 2);
                nextIndex = x*1.8 + 1;
            }
            oldYLogs[x-lowLimit] = lastLog;
        }
        
        for (int j = 0; j < maxIndex; j++){//for each prime
            const long prime = oldPrime[j];
            const long root = oldRoot[j];
            const auto primeLog = oldPrimeLogs[j];
            //cout<<j<<" "<<endl;
            auto i = lowLimit - lowLimit%prime + root;
            if (i<lowLimit) i+=prime;
            
            //auto i = oldi[j];
            //assert(i2==i);
            //cout<<i<<" "<<i2<<endl;
            for (; i<highLimit; i += prime){
                oldYLogs[i-lowLimit] -= primeLog;
            }
            //oldi[j] = i;
        }
        for (auto i = lowLimit; i<highLimit; i++){
            if (oldYLogs[i-lowLimit] < maxPrimeLog) {
                mpz_class y = multiPolynomial(i);
                
                vector<long> v;
                for(auto p = 0; p < primeBase.size(); p++) {
                    while(mpz_divisible_ui_p(y.get_mpz_t(), primeBase[p].first)) {
                        mpz_divexact_ui(y.get_mpz_t(), y.get_mpz_t(), primeBase[p].first);
                        //                            bitset ^= (((mpz_class)1) << p);
                        v.push_back(p);
                    }
                }
                if (y == 1) {
                    //Lock this
                    while (lock.test_and_set(std::memory_order_acquire));  // acquire lock

                    if (sieveCount >= bitWidth) {
                        lock.clear(std::memory_order_release);
                        break;
                    }
                    
                    for (auto p: v) {
                        flipBit(bitsets[p], sieveCount);
                    }
                    oldX.push_back(i);
                    oldY.push_back((quotientSqrt + i)*(quotientSqrt + i) - quotient);
                    sieveCount++;
                    lock.clear(std::memory_order_release);               // release lock

                    //end lock
                }
            }
            
        }
    }
    
    FactorNumber quadraticSieve() {


        if (quotient == 1) {
            return FactorNumber(this->number, this->factors, this->quotient);
        }
        
        if (mpz_perfect_power_p(quotient.get_mpz_t())) {
            mpz_class root = 2;
            mpz_class power = 1337;
            auto exponent = 2;
            while (true) {
                mpz_root(root.get_mpz_t(), quotient.get_mpz_t(), exponent);
                mpz_pow_ui(power.get_mpz_t(), root.get_mpz_t(), exponent);
                if (power == quotient) {
                    vector<pair<mpz_class, long>> factors(this->factors);
                    factors.push_back(pair<mpz_class, long>(root, exponent));
                    return FactorNumber(this->number, factors, 1);
                }
                exponent++;
            }
            
            
        }
        
        // generating prime base
        auto logger = Logger();
        logger.log("Generating prime base");
        auto primeBase = generatePrimeBase2();
        
        //        cout << "Primes" << endl;
        //        for (auto p: primeBase) {
        //            cout << p.first << endl;
        //        }
        //        auto primeBase2 = generatePrimeBase2();
        //
        //        cout<<"Sizes are "<<primeBase.size()<<" "<<primeBase2.size()<<endl;
        //
        //        for (auto p = 0; p < primeBase.size(); p++) {
        //            auto primePair = primeBase[p];
        //            cout << "Prime number: " << primePair.first << endl;
        //            printVector(primePair.second);
        //        }
        //        cout<<endl;
        //        for (auto p = 0; p < primeBase2.size(); p++) {
        //            auto primePair = primeBase2[p];
        //            cout << "Prime number: " << primePair.first << endl;
        //            printVector(primePair.second);
        //        }
        //        cout<<endl<<endl;
        
        //
        auto count = primeBase.size() * 10;
        logger.log("Generating Y");
        // Generating Y
        
        auto rows2 = primeBase.size();
        auto columns2 = count;
        
        logger.log("Transposing matrix " + to_string(rows2) + "x" + to_string(columns2));
        
        
        //        for (long x = 0; x < count; x++) {
        //            mpz_class q = (quotientSqrt + x)*(quotientSqrt + x) - quotient;
        //            y.push_back(q);
        //        }
        
        logger.log("Sieving");
        // Sieving
        vector<mpz_class> oldY;
        vector<mpz_class> oldX;
        
        long bitWidth = primeBase.size()*1.05+100;

        vector<bitarray> bitsets(primeBase.size(),bitWidth);
        bitarray tmpBitArray(bitWidth);
        
        vector<long> oldi;
        vector<long> oldPrime;
        vector<float> oldPrimeLogs;
        //vector<float> oldYLogs;
        vector<long> oldRoot;

        
        long index=0;
        for (long long p = 0; p < primeBase.size(); p++) {
            auto primePair = primeBase[p];
            for (long root: primePair.second) {
                oldi.push_back(root);
                oldPrimeLogs.push_back(log(primePair.first) / log(2));
                oldPrime.push_back(primePair.first);
                oldRoot.push_back(root);
                index++;
            }
        }
        long maxIndex = index;
        
        long sieveCount = 0;
        long lowLimit=0;
        int chunkSize = 1024*32*64;
        long highLimit = chunkSize;
        //const auto maxPrimeLog = oldPrimeLogs[oldPrimeLogs.size()-1];
        //float nextIndex = 0;
        //float lastLog = 0;
        std::atomic_flag lock = ATOMIC_FLAG_INIT;
        
        int threadCount = 8;
        thread * threads = new thread[threadCount];
        mpz_class debugCount = 0;

        
        for (int i=0;i<threadCount;i++){
            threads[i] = thread(
            [&](long oldYLogs[]){
                while (lock.test_and_set(std::memory_order_acquire));  // acquire lock
                while (sieveCount < bitWidth) {
                    if (debugCount++%100 == 0)
                        cout<<"New chunk "<<lowLimit<<" "<<sieveCount<<"/"<<primeBase.size()<<endl;
                    lock.clear(std::memory_order_release);               // release lock
                    sieve(lowLimit, highLimit, oldPrimeLogs, primeBase, oldPrime, oldRoot, maxIndex, bitWidth, oldYLogs, oldY, oldX, sieveCount, bitsets, lock, oldi);
                    while (lock.test_and_set(std::memory_order_acquire));  // acquire lock
                    lowLimit = highLimit;
                    highLimit += chunkSize;
                    
                }
                //delete [] oldYLogs;
                lock.clear(std::memory_order_release);               // release lock
                delete[] oldYLogs;
            },new long[highLimit-lowLimit]);
        }
        for (int i = 0; i < threadCount; i++){
            threads[i].join();
        }
        cout << "Sieve size " << sieveCount << endl;
        cout << "Prime base sise " << primeBase.size() << endl;
        
        auto rows = primeBase.size();
        auto columns = oldY.size();
        
        cout << "Matrix" << endl;
        
        //        cout << "Bit sets" << endl;
        //        printBitVector(bitsets, rows);
        //        cout << "Matrix" << endl;
        //        printBitVector(matrix, columns);
        
        //vector<bitarray> matrix = bitsets;
        mpz_class bitCount=0;
        //for (int i=0;i<rows;i++)
        //    bitCount += countOnes(bitsets[i]);
        //cout<<"bitbitbitbitbitbitbitbitbitbi "<<bitCount<<" bitbitbitibtibitbi"<<endl;
        logger.log("Gauss elmination");
        cout<<"rows: "<<rows<<" columns"<<columns<<endl;
        // Echelon matrix
        long i = 0;
        long j = columns-1;
        long lastNonZeroRow = -1;
        while (i < rows && j >= 0) {
            for (auto row = i+1; row < rows; row++) {
                if (isBitSet(bitsets[row], j)) {
                    swap(bitsets[i], bitsets[row]);
                    break;
                }
            }
            if (isBitSet(bitsets[i], j)) {
                lastNonZeroRow = i;
                for (auto row = i + 1; row < rows; row++) {
                    if (isBitSet(bitsets[row], j)) {
                        mpn_xor_n(bitsets[row].limb, bitsets[row].limb, bitsets[i].limb, bitsets[row].limbCount);
//                        matrix[row] ^= matrix[i];
                    }
                }
                i++;
            }
            j--;
        }
        
        logger.log("Calculating solution");
        auto iterations = 1000;
        while (iterations > 0) {
            iterations--;
            gmpRandom.get_z_bits(columns);
            bitarray solution(bitWidth);
            mpn_random(solution.limb,solution.limbCount);
            //            solution = ~solution;
            
            //            cout << "Solution start value "  << endl;
            //            printBits(solution, columns);
            
            long row = lastNonZeroRow;
            while (row >= 0) {
                mpn_and_n(tmpBitArray.limb, solution.limb, bitsets[row].limb, solution.limbCount);
                if (countOnes(tmpBitArray) % 2 == 1) {
                    flipBit(solution, leftMostOne(bitsets[row]));
                }
                //mpn_and_n(tmpBitArray.limb, solution.limb, matrix[row].limb, solution.limbCount);
                //assert(countOnes(tmpBitArray) % 2 == 0);

                row--;
            }
            
            //            cout << "Solution final" << endl;
            //            printBits(solution, columns);
            
//            for (auto row = 0; row < rows; row++) {
//                mpn_and_n(tmpBitArray.limb, solution.limb, matrix[row].limb, solution.limbCount);
//                if (countOnes(tmpBitArray) % 2 != 0) {
//                    //printBitVector(matrix, columns);
//                    cout << "LEFT MOST " << leftMostOne(matrix[0]) << endl;
//                    cout << "LEFT MOST " << leftMostOne(matrix[row]) << endl;
//                    cout << "Rows " << rows << endl;
//                    cout << "Row " << row << endl;
//                    cout << "Last non zero row " << lastNonZeroRow << endl;
//                    cout << "ERROR A'LA CARTÃ©" << endl;
//                }
//            }
            
            mpz_class Y = 1;
            mpz_class X = 1;
            for (auto i = 0; i < columns; i++) {
                if (isBitSet(solution, i)) {
                    Y *= oldY[i];
                    X *= (oldX[i] + quotientSqrt);
                }
            }
            if (msqrtceiling(Y) != msqrtfloor(Y)) {
                cout << "ERROR ultra bug in roots" << endl;
            }
            //            cout << "Y^2: " << Y << endl;
            
            Y = msqrtceiling(Y);
            
            //            cout << "Y^2: " << Y << endl;
            //            cout << "X^2: " << X << endl;
            cout << "N: " << quotient << endl;
            cout << "X-Y: " << gcd((X-Y), quotient) << endl;
            cout << "X+Y: " << gcd((X+Y), quotient / gcd((X-Y), quotient)) << endl;
            
            vector<pair<mpz_class, long>> factors(this->factors);
            
            mpz_class x = quotient;
            
            //        if (X % quotient == Y || X % quotient == -Y) {
            //            continue;
            //        }
            
            
            auto commonFactors = gcd(X-Y,quotient);
            if (commonFactors != 1  && commonFactors != quotient) {
                if (!isPrime(commonFactors) && false) {
                    cout << "Re seiving " << commonFactors << endl;
                    vector<pair<mpz_class, long>> v;
                    auto number = FactorNumber(commonFactors, v, commonFactors);
                    number = number.quadraticSieve();
                    cout << " New factors " << commonFactors <<  " " << number.factors.size() << endl;
                    //                    printVector(number.factors);
                    factors.insert(factors.end(), number.factors.begin(), number.factors.end());
                    x = number.quotient;
                } else {
                    factors.push_back(pair<mpz_class, long>(commonFactors, 1));
                    x /= (commonFactors);
                }
            }
            commonFactors = gcd(X+Y,x);
            if (commonFactors != 1  && commonFactors != quotient) {
                if (!isPrime(commonFactors) && false) {
                    cout << "Re seiving " << commonFactors << endl;
                    vector<pair<mpz_class, long>> v;
                    auto number = FactorNumber(commonFactors, v, commonFactors);
                    number = number.quadraticSieve();
                    factors.insert(factors.end(), number.factors.begin(), number.factors.end());
                    x = number.quotient;
                } else {
                    factors.push_back(pair<mpz_class, long>(commonFactors, 1));
                    x /= (commonFactors);
                }
            }
            
            if (factors.size() > this->factors.size()) {
                logger.log("Solution found");
                return FactorNumber(number, factors, x);
            } else {
                logger.log("No solution found " + to_string(iterations));
            }
        }
        return FactorNumber(number, factors, quotient);
    }
    
    vector<pair<long, vector<long>>> generatePrimeBase(vector<long> primes) {
        vector<pair<long, vector<long>>> primeBase;
        for (auto prime: primes) {
            if (prime > B) {
                break;
            }
            auto mod = quotient % prime;
            for (long i = 0; i < prime; i++) {
                if ((i*i) % prime == mod) {
                    if (primeBase.empty() || primeBase.back().first != prime) {
                        pair<long, vector<long>> primeRoots(prime, vector<long>());
                        primeBase.push_back(primeRoots);
                    }
                    mpz_class tmp = (((i - quotientSqrt) % prime + prime) % prime);
                    primeBase.back().second.push_back(tmp.get_si());
                }
            }
        }
        return primeBase;
    }
    
    
    
    vector<pair<long, vector<long>>> generatePrimeBase2() {
        vector<pair<long, vector<long>>> primeBase;
        
        for (auto prime: primes) {
            if (prime > B) {
                break;
            }
            //            cout << "legender " << prime << " " << legendreSymbol(quotient, prime) << endl;
            if (prime == 2 || legendreSymbol(quotient, prime) == 0) {
                vector<long> vec2{prime};
                auto base = generatePrimeBase(vec2);
                primeBase.push_back(base.back());
                continue;
            }
            if (prime<3) continue;
            if (legendreSymbol(quotient, prime) != 1) continue;
            
            
            mpz_class s = 0;
            mpz_class q = prime-1;
            
            while ( q%2 == 0){
                q /= 2;
                s++;
            }
            
            if (s==1){
                mpz_class a = powMod(quotient, (prime+1)/4, prime);
                mpz_class b = (-a+prime) % prime;
                mpz_class root1 = (((a - quotientSqrt) % prime + prime) % prime);
                mpz_class root2 = (((-a - quotientSqrt) % prime + prime) % prime);
                primeBase.push_back(pair<int, vector<long>>(prime,vector<long>{root1.get_si(),root2.get_si()}));
                continue;
            }
            
            //Hitta z
            mpz_class z = 2;
            while (legendreSymbol(z, prime) != -1)
                z++;
            
            mpz_class c = powMod(z, q, prime);
            mpz_class r = powMod(quotient, (q+1)/2, prime);
            mpz_class t = powMod(quotient, q, prime);
            mpz_class m = s;
            
            while (t % prime != 1){
                int i = 1;
                while (powMod(t, powNoMod(2, i), prime) != 1)
                    i++;
                auto b = powMod(c, powNoMod(2, m.get_si()-i-1), prime);
                r = (r * b) % prime;
                t = (t * b * b) % prime;
                c = (b * b) % prime;
                m = i;
            }
            //Solution is now p+r and p-r
            
            //            cout << "roots for 11 " << r << " " << prime -r << " " << prime + r;
            mpz_class root1 = (((r - quotientSqrt) % prime + prime) % prime);
            mpz_class root2 = ((((prime-r) - quotientSqrt) % prime + prime) % prime);
            
            primeBase.push_back(pair<int, vector<long>>(prime,vector<long>{root1.get_si(),root2.get_si()}));
        }
        
        return primeBase;
    }
};
ofstream primeFile;
bool factorize(mpz_class n,int id) {
    vector<pair<mpz_class, long>> v;
//    cout<<"doing number "<<n<<endl;
    auto number = FactorNumber(n, v, n);
    number = number.primalDivision(1000).quadraticSieve();
    vector<pair<mpz_class, long>> primeFactors;
    vector<pair<mpz_class, long>> factors(number.factors);
    if (number.quotient != 1) {
        pair<mpz_class, long> p(number.quotient, 1);
        factors.push_back(p);
    }
    for (int i = 0; i < factors.size(); i++) {
        cout << "Workin on factor " << factors[i].first << endl;
        if (!isPrime(factors[i].first)) {
            vector<pair<mpz_class, long>> v;
            auto number = FactorNumber(factors[i].first, v, factors[i].first).quadraticSieve();
            for (auto factor: number.factors) {
                pair<mpz_class, long> v(factor.first, factor.second * factors[i].second);
                factors.push_back(v);
            }
            if (number.quotient != 1) {
                pair<mpz_class, long> p(number.quotient, 1);
                factors.push_back(p);
            }
        } else {
            primeFactors.push_back(factors[i]);
        }
    }
    cout << number.number << " factorized "<<id << endl;
    mpz_class sum = 1;
    bool isPrimes = true;
    primeFile << number.number << endl;
    for (auto primeFactor: primeFactors) {
        isPrimes = isPrimes & isPrime(primeFactor.first);
        sum *= powNoMod(primeFactor.first, primeFactor.second);
        primeFile << primeFactor.first << "^" << primeFactor.second << "*";
        cout << primeFactor.first << "^" << primeFactor.second << "*";
    }
    assert(sum==n);
    assert(isPrimes);
    primeFile << endl;
    cout << endl;
    
    return sum == n && isPrimes;
}




int maxNumber=50;
int parts = 1;
int readIndex=0;
vector<pair<int,mpz_class>> numbers;

void threadStart(){
    cout<<"hej"<<endl;
    Logger l2;
    while(true){
        if (!(readIndex<numbers.size())) break;
        auto number = numbers[readIndex];
        readIndex++;
        l2.log("------STARTING FACTORIZING-----");
        //cout<<"Working on "<<number<<endl;
        factorize(mpz_class(number.second),number.first);
        l2.log("-----DONE FACTORIZING-------");
    }
    
}



int main(int argc, const char * argv[]) {
    
    primeFile.open("allFile");
    
    mpz_class n("9108020935");
    mpz_class big;
    mpz_pow_ui(big.get_mpz_t(), ((mpz_class)10).get_mpz_t(), 45);
    n*=big;
    
    //for (int i=1;i<=1;i++)
    //    numbers.push_back(pair<int, mpz_class>(i, n+i));
    

    //auto paulCompleted ={1, 82, 80, 64, 63, 58, 57, 56, 53, 51, 50, 49, 48, 47, 44, 43, 38, 35, 34, 37, 32, 30, 29, 28, 25, 24, 21, 22, 23, 20, 19, 17, 14, 13, 12, 10, 9, 6, 4, 2, 75, 84, 55, 36, 79, 88, 11, 45, 91, 85, 27, 11, 87, 95, 94, 97};
    //for (int i=1;i<=100;i++)
    //    if (find(paulCompleted.begin(), paulCompleted.end(), i) == paulCompleted.end())
    //        numbers.push_back(pair<int, mpz_class>(i, n+i));
    
    auto todo = {5};//,7,15,16,26,31,33,40,41,52,54,59,60,61,62,65,66,68,70,72,73,77,78,81,83,86,89,90,92,96,100};
    for (int i=1;i<=100;i++)
        if (find(todo.begin(), todo.end(), i) != todo.end())
            numbers.push_back(pair<int, mpz_class>(i, n+i));
    
    
//    
//    n = "9011221992";
//    n *= big;
//    
//    for (int i=1;i<=100;i++)
//        numbers.push_back(pair<int, mpz_class>(i, n+i));

    
    thread **threads = (thread**)malloc(sizeof(thread)*parts);
    
    for (int c=0; c<parts; c++){
        threads[c] = new thread(threadStart);
        chrono::milliseconds duration(100);
        this_thread::sleep_for(duration);
        
        
    }
    for (int c=0; c<parts; c++){
        threads[c]->join();
    }
    
    return 0;
}

