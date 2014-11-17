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
using namespace std;

void setBit(mpz_class &x, long i) {
    x |= (((mpz_class)1) << i);
}

void unsetBit(mpz_class &x, long i) {
    x &=~ (((mpz_class)1) << i);
}

bool isBitSet(mpz_class x, long i) {
    return ((x >> i) & 1) == 1 ? true : false;
}


template <typename T>
void printVector(std::vector<T> vector) {
    for (auto t: vector) {
        std::cerr << t << " ";
    }
    std::cerr << std::endl;
}

void printBitVector(std::vector<mpz_class> vector, size_t length) {
    for (auto v: vector) {
        for (int i = 0; i < length; i++) {
            cout << isBitSet(v, i);
        }
        cout << endl;
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



void writePrimesToFile(string fileName, vector<long> primes) {
    ofstream file;
    file.open(fileName);
    for (auto prime: primes) {
        file << prime << endl;
    }
}

bool isPerfectPower(mpz_class number) {
    return mpz_perfect_power_p(number.get_mpz_t()) == 1;
}

mpz_class legendreSymbol(mpz_class n,mpz_class p){
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
    writePrimesToFile(primeFileName, primes);
    return primes;
}

auto primes = generatePrimes(1e8);

mpz_class pollard(mpz_class n, long startValue, mpz_class limit) {
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

long a = 1;

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

bool hasRoot(mpz_class x, long p){
    mpz_class tmp;
    mpz_powm_ui(tmp.get_mpz_t(),x.get_mpz_t(),(p-1)/2,mpz_class(p).get_mpz_t());
    //        mpz_pow_ui(tmp.get_mpz_t(),x.get_mpz_t(),(p-1)/2);
    //        mpz_mod_ui(tmp.get_mpz_t(),tmp.get_mpz_t(),p);
    return tmp == 1;
}


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
        double n = number.get_d();
        B = 3*exp(0.5*sqrt(log(n)*log(log(n)))) + 10;
        //        print();
        
        if (isPrime(quotient)) {
            cout << "Quotient is prime" << endl;
            this->factors.push_back(pair<mpz_class, long>(quotient, 1));
            this->quotient = 1;
        }
    }
    
    FactorNumber pollardish(mpz_class limit) {
        cout << "Pollardish" << endl;
        vector<pair<mpz_class, long>> factors(this->factors);
        auto quotient = this->quotient;
        for (auto startValue = 2; startValue < 5; startValue++) {
            if (isPrime(quotient)) {
                cout << "IS PRIME" << endl;
                break;
            }
            auto factor = pollard(quotient, startValue, limit);
            if (factor != NULL) {
                cout << "Pollard found " << factor << endl;
                factors.push_back(std::pair<mpz_class, long>(factor, 1));
                quotient /= factor;
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
    
    
    FactorNumber primalDivision() {
        cout << "Primal division" << endl;
        mpz_class x = quotient;
        vector<pair<mpz_class, long>> factors;
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
    
    mpz_class mpow(mpz_class x, long a) {
        mpz_class pow;
        mpz_pow_ui(pow.get_mpz_t(), x.get_mpz_t(), a);
        return pow;
    }
    
    FactorNumber quadraticSieve() {
        // generating prime base
        
        auto primeBase = generatePrimeBase();
        auto primeBase2 = generatePrimeBase2();
        
        cout<<"Sizes are "<<primeBase.size()<<" "<<primeBase2.size()<<endl;
        
        for (auto p = 0; p < primeBase.size(); p++) {
            auto primePair = primeBase[p];
            cout << "Prime number: " << primePair.first << endl;
            printVector(primePair.second);
        }
        cout<<endl;
        for (auto p = 0; p < primeBase2.size(); p++) {
            auto primePair = primeBase2[p];
            cout << "Prime number: " << primePair.first << endl;
            printVector(primePair.second);
        }
        cout<<endl<<endl;
        auto count = 10 * primeBase.size();
        
        // Generating Y
        vector<mpz_class> y;
        for (long x = 0; x < count; x++) {
            mpz_class q = (quotientSqrt + x)*(quotientSqrt + x) - quotient;
            y.push_back(q);
        }
        
        // Sieving
        vector<mpz_class> oldY(y);
        vector<mpz_class> bitsets(y.size(), 0);
        vector<vector<long>> bitsetsCount(y.size(), vector<long>(primeBase.size(), 0));
        
        for (long long p = 0; p < primeBase.size(); p++) {
            auto primePair = primeBase[p];
            for (long root: primePair.second) {
                for (long i = root; i < y.size(); i += primePair.first) {
//                    cout << y[i] << " " << primePair.first << " " << y[i] % primePair.first << endl;
                    while (y[i] % primePair.first == 0) {
                        y[i] /= primePair.first;
                        bitsets[i] ^= (((mpz_class)1) << p);
                        bitsetsCount[i][p]++;
//                        if (bitsetsCount[i][p] % 2 != ((bitsets[i] >> p) & 1)) {
//                            cout << "ULTRA ERROR BITWISE STUFF DOES NOT WORK" << endl;
//                        }
                    }
                }
            }
        }
        
        cout << "Bitsets" << endl;
        printVector(bitsets);
        
        cout << "Y" << endl;
        printVector(y);
        
        // Finding Y's that factored in our base
        for (long i = 0; i < y.size(); i++) {
            if (y[i] == 1) {
                for (long p = 0; p < primeBase.size(); p++) {
                    cout << ((bitsets[i] >> p) & 1) << " ";
                }
                cout << endl;
                cout << "Found something good " << i << endl;
            } else {
                bitsets[i] = 0;
            }
        }
        

        
        auto rows = primeBase.size();
        auto columns = y.size();
        vector<mpz_class> matrix(rows, mpz_class(0));
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < columns; column++) {
                if (isBitSet(bitsets[column], row)) {
                    setBit(matrix[row], column);
                }
            }
        }
        cout << "Bitsets" << endl;
//        printBitVector(bitsets, rows);
        cout << "Transposed?" << endl;
        printBitVector(matrix, columns);
        // Ugly finding solution
        
        // Echelon matrix
        int i = 0;
        int j = 0;
        while (i < rows && j < columns) {
            for (int row = i+1; row < rows; row++) {
                if (isBitSet(matrix[row], j)) {
                    swap(matrix[i], matrix[row]);
                    break;
                }
            }
            if (isBitSet(matrix[i], j)) {
                for (int row = i + 1; row < rows; row++) {
                    if (isBitSet(matrix[row], j)) {
                        matrix[row] ^= matrix[i];
                    }
                }
                i++;
            }
            j++;
        }
        
        
        cout << "Echelon?" << endl;
        printBitVector(matrix, columns);
        
//        int a = 0;
//        int b = 0;
//        while(a % quotient == b % quotient || a % quotient == (- b) % quotient + quotient) {
//            
//        }
        
        
        for (long i = 0; i < bitsets.size(); i++) {
            for (long j = i+1; j < bitsets.size(); j++) {
                for (long k = j+1; k < bitsets.size(); k++) {
                    if (bitsets[i] != 0 && bitsets[j] != 0 && bitsets[k] != 0) {
                        auto zero = bitsets[i]^bitsets[j]^bitsets[k];
                        if (zero == 0) {
                            cout << "Perfect match" << endl;
                            cout << oldY[i] << " " << oldY[j] << " " << oldY[k] << endl;
                            cout << bitsets[i] << " " << bitsets[j] << " " << bitsets[k] << endl;
//                            printVector(bitsetsCount[i]);
//                            printVector(bitsetsCount[j]);
//                            printVector(bitsetsCount[k]);
//                            cout << bitsetsCount[i][ << " " << bitsetsCount[j] << " " << bitsetsCount[k] << endl;
                            
                            for (auto p: primeBase) {
                                cout << p.first << " ";
                            }
                            cout << endl;
//                            cout << primeBase[i].first << " " << primeBase[j].first << " " << primeBase[k].first << endl;
                            mpz_class Y = 1;
                            
                            for (long index: {i,j,k}) {
                                Y *= oldY[index];
                            }
                            if (msqrtceiling(Y) != msqrtfloor(Y)) {
                                cout << "ERROR ultra bug in roots";
                            }
                            Y = msqrtceiling(Y);
                            
                            cout << bitsets[i] << " " << bitsets[j] << " " << bitsets[k] << endl;
                            cout << "Y^2: " << Y << endl;
                            cout << i << " " << j << " " << k << endl;
                            cout << i << " " << j << " " << k << endl;
                            cout << i + quotientSqrt << " " << j + quotientSqrt << " " << k + quotientSqrt << endl;
                            mpz_class X = (i + quotientSqrt)*(j + quotientSqrt)*(k + quotientSqrt);
                            cout << "X^2: " << X << endl;
                            cout << "N: " << quotient << endl;
                            cout << "X-Y: " << gcd((X-Y), quotient) << endl;
                            cout << "X+Y: " << gcd((X+Y), quotient) << endl;
                            
                            vector<pair<mpz_class, long>> factors(this->factors);
                            
                            mpz_class x = quotient;
                            
                            if (X % quotient == Y || X % quotient == -Y) {
                                continue;
                            }
                            
                            auto commonFactors = gcd(X-Y,quotient);
                            if (commonFactors != 1  && commonFactors != quotient) {
                                factors.push_back(pair<mpz_class, long>(commonFactors, 1));
                                x /= (commonFactors);
                            }
                            commonFactors = gcd(X+Y,x);
                            if (commonFactors != 1  && commonFactors != quotient) {                                factors.push_back(pair<mpz_class, long>(commonFactors, 1));
                                x /= (commonFactors);
                            }
                            
                            if (factors.size() > this->factors.size()) {
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
    
    vector<pair<long, vector<long>>> generatePrimeBase() {
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
            if (prime<3) continue;
            if (legendreSymbol(quotient, prime) != 1) continue;

            
            mpz_class s = 0;
            mpz_class q = prime-1;
            
            while ( q%2 == 0){
                q /= 2;
                s++;
            }
            
            //Fast method but only works if s=1
            if (s==1){
                mpz_class a = powMod(quotient, (prime+1)/4, prime);
                mpz_class b = (-a+prime) % prime;
                primeBase.push_back(pair<int, vector<long>>(prime,vector<long>{a.get_si(),b.get_si()}));
                continue;
            }
            
            //Hitta z
            mpz_class z = 1;
            while (legendreSymbol(z, prime) != -1)
                z++;
            
            mpz_class c = powMod(z, q, prime);
            mpz_class r = powMod(quotient, (q+1)/2, prime);
            mpz_class t = powMod(quotient, q, prime);
            mpz_class m = s;
            
            while (t != 1 && t != -prime+1){
                auto t2 =powMod(t, 2, prime);
                int i = 1;
                while (powMod(t2, i, prime) != 1)
                    i++;
                
                auto b = powMod(c, powMod(2, m-i-1, prime), prime);
                
                r = (r * b) % prime;
                t = (t * b * b) % prime;
                c = (b * b) % prime;
                m = i;
            }
            //Solution is now p+r and p-r
            
            primeBase.push_back(pair<int, vector<long>>(prime,vector<long>{((mpz_class)(prime+r)).get_si(),((mpz_class)(prime-r)).get_si()}));
            
        }
        return primeBase;
    }
};

int main(int argc, const char * argv[]) {
    mpz_class n("9011221992");
    mpz_class big;
    mpz_pow_ui(big.get_mpz_t(), ((mpz_class)10).get_mpz_t(), 20);
    n *= big;
    n += 1;
    n = 12;
//    n = 15348;
//    n = 9011221992;
    vector<pair<mpz_class, long>> v;
    auto number = FactorNumber(n, v, n).quadraticSieve();
    number.internalCheck();
    number.print();
    return 0;
}