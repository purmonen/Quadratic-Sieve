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

void setBit(mpz_class &x, long i) {
    x |= (((mpz_class)1) << i);
}

void unsetBit(mpz_class &x, long i) {
    x &=~ (((mpz_class)1) << i);
}

void flipBit(mpz_class &x, long i) {
    x ^= (((mpz_class)1) << i);
}

bool isBitSet(mpz_class x, long i) {
    return ((x >> i) & 1) == 1 ? true : false;
}

long countOnes(mpz_class x) {
    mpz_class count = 0;
    while (x != 0) {
        count += x & 1;
        x >>= 1;
    }
    return count.get_si();
}

long leftMostOne(mpz_class x) {
    mpz_class count = -1;
    while (x != 0) {
        count++;
        x = x >> 1;
    }
    
    return count.get_si();
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
    cout << std::endl;
}

void printBits(mpz_class v, size_t length) {
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
//    writePrimesToFile(primeFileName, primes);
    return primes;
}

auto primes = generatePrimes(1e7);


mpz_class pollardLastX = 0;
mpz_class pollardLastY = 0;
mpz_class pollardLastStartValue = 0;
mpz_class pollard(mpz_class n, long startValue, mpz_class limit) {
    mpz_class x = startValue, y = startValue, d = 1;
    auto startTime = chrono::system_clock::now();
    int iterations = 0;
    while (d == 1) {
        iterations++;
        if (iterations > 10000) {
            iterations = 0;
            if (((chrono::system_clock::now() - startTime).count() / 1e6) > limit) {
                cout << "Hit pollard limit " << startValue << endl;
                break;
            }
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
        double n = number.get_d();
        B = 3*exp(0.5*sqrt(log(n)*log(log(n)))) + 30;
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
        
        for (auto startValue = 2; startValue <= 3; startValue++) {
            if (isPrime(quotient)) {
                cout << "IS PRIME" << endl;
                break;
            }
            auto factor = pollard(quotient, startValue, 10);
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
        if (quotient == 1) {
            return FactorNumber(this->number, this->factors, this->quotient);
        }

        // generating prime base
        

        
        

        
        auto logger = Logger();
        logger.log("Generating prime base");
        auto primeBase = generatePrimeBase(primes);
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
        auto count = primeBase.size() * 2;
        
        
        logger.log("Generating Y");
        // Generating Y
        
        auto rows2 = primeBase.size();
        auto columns2 = count;
        
        logger.log("Transposing matrix " + to_string(rows2) + "x" + to_string(columns2));
        
        vector<mpz_class> y;
        for (long x = 0; x < count; x++) {
            mpz_class q = (quotientSqrt + x)*(quotientSqrt + x) - quotient;
            y.push_back(q);
        }
        
        logger.log("Sieving");
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
        
        // Finding Y's that factored in our base
        for (long i = 0; i < y.size(); i++) {
            if (y[i] == 1) {
//                for (long p = 0; p < primeBase.size(); p++) {
//                    cout << ((bitsets[i] >> p) & 1) << " ";
//                }
//                cout << endl;
//                cout << "Found something good " << i << endl;
            } else {
                bitsets[i] = 0;
            }
        }

        auto rows = primeBase.size();
        auto columns = y.size();
        logger.log("Transposing matrix " + to_string(rows) + "x" + to_string(columns));
        vector<mpz_class> matrix(rows, mpz_class(0));
        for (auto row = 0; row < rows; row++) {
            for (auto column = 0; column < columns; column++) {
                if (isBitSet(bitsets[column], row)) {
                    setBit(matrix[row], column);
                }
            }
        }
        
        
        
//        cout << "Bit sets" << endl;
//        printBitVector(bitsets, rows);
//        cout << "Matrix" << endl;
//        printBitVector(matrix, columns);
        
        
        logger.log("Gauss elmination");
        // Echelon matrix
        long i = 0;
        long j = columns-1;
        long lastNonZeroRow = -1;
        while (i < rows && j >= 0) {
            for (auto row = i+1; row < rows; row++) {
                if (isBitSet(matrix[row], j)) {
                    swap(matrix[i], matrix[row]);
                    break;
                }
            }
            if (isBitSet(matrix[i], j)) {
                lastNonZeroRow = i;
                for (auto row = i + 1; row < rows; row++) {
                    if (isBitSet(matrix[row], j)) {
                        matrix[row] ^= matrix[i];
                    }
                }
                i++;
            }
            j--;
        }
        

        
        
        logger.log("Calculating solution");
        vector<mpz_class> solutionMatrix(matrix);
        gmpRandom.get_z_bits(columns);
        mpz_class solution = 0;//gmpRandom.get_z_bits(columns);
        solution = ~solution;

        for (int i = 0; i < columns; i++) {
            if (y[i] != 1) {
                unsetBit(solution, i);
            }
        }
        
        long row = rows-1;
        while (row >= 0) {
            if (countOnes(solution & matrix[row]) % 2 == 1) {
                flipBit(solution, leftMostOne(matrix[row]));
            }
            if (countOnes(solution & matrix[row]) % 2 == 1) {
                cout << "WHAT?!?!?!?!? SHOULD BE 0" << endl;
            }
            row--;
        }
        
        for (auto row = 0; row < rows; row++) {
            if (countOnes(solution & matrix[row]) % 2 != 0) {
                printBitVector(matrix, columns);
                cout << "LEFT MOST " << leftMostOne(matrix[0]) << endl;
                cout << "LEFT MOST " << leftMostOne(matrix[row]) << endl;
                cout << "Rows " << rows << endl;
                cout << "Row " << row << endl;
                cout << "Last non zero row " << lastNonZeroRow << endl;
                cout << "ERROR A'LA CARTé" << endl;
            }
        }
        
        mpz_class Y = 1;
        mpz_class X = 1;
        logger.log("Multiplying large numbers");
        for (auto i = 0; i < columns; i++) {
            if (isBitSet(solution, i)) {
                Y *= oldY[i];
                X *= ((i) + quotientSqrt);
            }
        }
        if (msqrtceiling(Y) != msqrtfloor(Y)) {
            cout << "ERROR ultra bug in roots" << endl;
        }
        cout << "Y^2: " << Y << endl;

        Y = msqrtceiling(Y);

        cout << "Y^2: " << Y << endl;
        cout << "X^2: " << X << endl;
        cout << "N: " << quotient << endl;
        cout << "X-Y: " << gcd((X-Y), quotient) << endl;
        cout << "X+Y: " << gcd((X+Y), quotient) << endl;
        
        vector<pair<mpz_class, long>> factors(this->factors);
        
        mpz_class x = quotient;
        
//        if (X % quotient == Y || X % quotient == -Y) {
//            continue;
//        }
        
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
            logger.log("Solution found");
            return FactorNumber(number, factors, x);
        } else {
            logger.log("No solution found");
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
            cout << "legender " << prime << " " << legendreSymbol(quotient, prime) << endl;
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

bool factorize(mpz_class n) {
    vector<pair<mpz_class, long>> v;
    auto number = FactorNumber(n, v, n);
    number = number.quadraticSieve();
    number.internalCheck();
    number.print();
    return number.quotient == 1;
}

int main(int argc, const char * argv[]) {

    mpz_class n("9011221992");
    mpz_class big;
    mpz_pow_ui(big.get_mpz_t(), ((mpz_class)10).get_mpz_t(), 10);
    n *= big;
    n += 0;
    n = 15350;
    int count = 0;

//    for (int i = 1; i < 100; i++) {
        if (factorize(n)) {
            count++;
        }
//    }
    cout << "Factored " << count << endl;
    return 0;
}