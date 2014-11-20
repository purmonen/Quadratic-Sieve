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

void setBit(mpz_class &x, long i) {
    x |= (((mpz_class)1) << i);
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


#define GMP_LIMB_BITS_SHIFTER 6
struct bitarray{
    mp_limb_t* limb;
    mp_size_t limbCount;
    int id=rand();
    
    bitarray(long bitCount){
        assert(bitCount!=0);
        cout<<"creating with normal constructor"<<endl;
        limbCount = (bitCount>>GMP_LIMB_BITS_SHIFTER)+1;
        cout<<"ID"<<id<<"allocating "<<limbCount<<" for "<<bitCount<<"bits"<<endl;
        limb = new mp_limb_t[limbCount];
        memset(limb, 0, limbCount*GMP_LIMB_BITS>>3);
    }
    
    bitarray(bitarray const &b){
        cout<<"creating with copy constructor"<<endl;
        limbCount = b.limbCount;
        limb = new mp_limb_t[limbCount];
        memcpy(limb, b.limb, limbCount*8);
    }
    
    ~bitarray(){
        cout<<"Destroy "<<this<<" delete "<<limb<<endl;
        delete[] limb;
    }
};

void swap(bitarray &a, bitarray &b){
    swap(a.limb,b.limb);
}

void setBit(bitarray &x, const long &i) {
    const long limb = i >> GMP_LIMB_BITS_SHIFTER;
    const long bitInLimb = i & 0b11111111ull;
    x.limb[limb] |= 1ull << bitInLimb;
}

void unsetBit(bitarray &x, const long &i) {
    const long limb = i >> GMP_LIMB_BITS_SHIFTER;
    const long bitInLimb = i & 0b11111111ull;
    x.limb[limb] &= ~(1ull << bitInLimb);
}

void flipBit(bitarray &x, const long &i) {
    const long limb = i >> GMP_LIMB_BITS_SHIFTER;
    const long bitInLimb = i & 0b11111111ull;
    x.limb[limb] ^= 1ull << bitInLimb;
}

inline const bool isBitSet(const bitarray &x, const long &i) {
    const long limb = i >> GMP_LIMB_BITS_SHIFTER;
    const long bitInLimb = i & 0b11111111ull;
    return (x.limb[limb]>>bitInLimb) & 1;
}

long countOnes(bitarray &x) {
    return mpn_popcount(x.limb,x.limbCount);
}
long leftMostOne(bitarray &bitset) {
    
    auto i = bitset.limbCount-1;
    for (; i>=0; i--) {
        if (bitset.limb[i] != 0) break;
    }
    cout<<"found in limb "<<i<<endl;
    auto x = bitset.limb[i];
    long count = -1;
    while (x != 0) {
        count++;
        x = x >> 1;
    }
    
    return (i<<GMP_LIMB_BITS_SHIFTER) + count;
}





int main(int argc, const char * argv[]) {
    bitarray arr1(100);
    for (int i=0;i<100;i++)
        setBit(arr1, i);
    
    bitarray arr2(100);
    for (int i=0;i<100;i+=2)
        setBit(arr2, i);

    swap(arr1,arr2);
    
    cout<<endl;
    for (int i=0;i<128;i++)
        cout<<isBitSet(arr1, i);
    cout<<endl;
    for (int i=0;i<128;i++)
        cout<<isBitSet(arr2, i);
    cout<<endl;
    
//    mpz_class arr2;
//    srand(5);
//    int bit =rand()%10;
//    setBit(arr, bit);
//    setBit(arr2,bit);
//    bit =rand()%10;
//    setBit(arr, bit);
//    setBit(arr2,bit);
//
//    cout<<leftMostOne(arr)<<endl;
//    cout<<leftMostOne(arr2)<<endl;
    
    return 1;
}






