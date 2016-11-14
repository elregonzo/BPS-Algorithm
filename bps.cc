#include <iostream>
#include <vector>
#include <cmath>
#define TWEAK_BASE 2LL
#define TWEAK_LEFT_LENGTH 32LL

using namespace std;


long long square_and_multiply(long long a, long long b){
    long long result = 1;

    while ( b > 0 ){
        if ( b % 2 == 1 )
            result *= a;
        b /= 2;
        a *= a;
    }

    return result;
}

/*
    Function that represents the Internal Block Cipher for the BPS cryptosystem, the parameters are described as follows:
    s: represents the cardinality of the alphabet
    b: represents the block length
    w: the number of rounds, by recommendation of the author should be 8 rounds
    X: The message to be cipher
    K: The key to be used in the F function
    T: Tweak to be used in the F function
*/
vector<long long> InternalBlockCipher(/*TODO: research how to insert a function here*/long long s, long long b
                        , long long w, vector<long long> X, vector<long long> K, unsigned long long T){

    long long Tr = T % square_and_multiply( TWEAK_BASE , TWEAK_LEFT_LENGTH);

    long long Tl = (T - Tr) / square_and_multiply( TWEAK_BASE , TWEAK_LEFT_LENGTH);

    long long l = ceil(b/2.0), r = floor( b/2.0 );

    vector< long long > L( w + 1 ,0) , R( w + 1 ,0), Y(b) ;
    for (int i = 0 ; i < l ; i++ )
        L[0] += X[i] * square_and_multiply( s , i );
    for (int i = 0 ; i < r ; i++ )
        R[0] += X[i+l] * square_and_multiply( s , i );

    for (int i = 0; i <= w - 1 ; i++ ){
        if ( i % 2 == 0 ) {
            L[i+1] = ( L[i] + FK((Tr ^ i)* square_and_multiply( 2 , f-32 ) + R[i]) ) % square_and_multiply( s , l );
            R[i+1] = R[i];
        }
        else{
            R[i+1] = ( R[i] + FK((Tl ^ i)* square_and_multiply( 2 , f-32 ) + L[i]) ) % square_and_multiply( s , r );
            L[i+1] = L[i];
        }

    }
    for (int i = 0 ; i <= l - 1 ; i ++  ){
        Y[i] = L[w] % s;
        L[w] = ( L[w] - Y[i] ) / s;
    }
    for (int i = 0 ; i <= r - 1 ; i ++  ){
        Y[ i + l ] = R[w] % s;
        R[w] = ( R[w] - Y[ i + l ] ) / s;
    }
    return Y ;
}

/*
    Function that represents the Internal Block Decipher for the BPS cryptosystem, the parameters are described as follows:
    s: represents the cardinality of the alphabet
    b: represents the block length
    w: the number of rounds, by recommendation of the author should be 8 rounds
    X: The message to be cipher
    K: The key to be used in the F function
    T: Tweak to be used in the F function
*/
vector<long long> InternalBlockCipher(/*TODO: research how to insert a function here*/long long s, long long b
                        , long long w, vector<long long> X, vector<long long> K, unsigned long long T){

    long long Tr = T % square_and_multiply( TWEAK_BASE , TWEAK_LEFT_LENGTH);

    long long Tl = (T - Tr) / square_and_multiply( TWEAK_BASE , TWEAK_LEFT_LENGTH);

    long long l = ceil(b/2.0), r = floor( b/2.0 );

    vector< long long > L( w + 1 ,0) , R( w + 1 ,0), Y(b) ;
    for (int i = 0 ; i < l ; i++ )
        L[w] += X[i] * square_and_multiply( s , i );
    for (int i = 0 ; i < r ; i++ )
        R[w] += X[i+l] * square_and_multiply( s , i );

    for (int i = w-1; i >= 0 ; i-- ){
        if ( i % 2 == 0 ) {
            L[i] = ( L[i+1] - FK((Tr ^ i)* square_and_multiply( 2 , f-32 ) + R[i+1]) ) % square_and_multiply( s , l );
            R[i] = R[i+1];
        }
        else{
            R[i] = ( R[i+1] - FK((Tl ^ i)* square_and_multiply( 2 , f-32 ) + L[i+1]) ) % square_and_multiply( s , r );
            L[i] = L[i+1];
        }

    }
    for (int i = 0 ; i <= l - 1 ; i ++  ){
        Y[i] = L[w] % s;
        L[w] = ( L[w] - Y[i] ) / s;
    }
    for (int i = 0 ; i <= r - 1 ; i ++  ){
        Y[ i + l ] = R[w] % s;
        R[w] = ( R[w] - Y[ i + l ] ) / s;
    }
    return Y ;

}

int main(){

    return 0;
}
