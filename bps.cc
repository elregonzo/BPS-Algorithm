#include <iostream>
#include <stdarg.h>
#include <obstack.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpxx.h>
#include <vector>

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
vector<long long> InternalBlockCipher(/*TODO: research how to insert a function here*/int s, int b
                        , int w, vector<int> X, vector<int> K, unsigned long long T){

    long long Tr = T % square_and_multiply( TWEAK_BASE , TWEAK_LEFT_LENGTH );

    long long Tl = (T - Tr) / square_and_multiply( TWEAK_BASE , TWEAK_LEFT_LENGTH);

    long long l = ceil(b/2.0), r = floor( b/2.0 );

    //Declare the BigIntegers to be used in each round 
    vector< mpz_t  > L( w + 1 ) , R( w + 1 ) ;

    //Initialize the BigIntegers nums
    for (int i = 0 ; i < w + 1 ; i++ )
    	mpz_inits( L[i] , R[i] , NULL );

    vector<long long> Y(b);

    //Transforms left branch from s base to decimal base
    for (int i = 0 ; i < l ; i++ ){
    	mpz_t exponentation;
    	mpz_init(exponentation); 
    	mpz_ui_pow_ui (exponentation, (unsigned long int) s , (unsigned long int) i );
    	mpz_addmul_ui (L[0], exponentation, (unsigned long int) X[i]);
    	mpz_clear(exponentation);
    }
    //Transforms right branch from s-base to decimal base
    for (int i = 0 ; i < r ; i++ ){
    	mpz_t exponentation;
    	mpz_init(exponentation); 
    	mpz_ui_pow_ui (exponentation, (unsigned long int) s , (unsigned long int) i );
    	mpz_addmul_ui (R[0], exponentation, (unsigned long int) X[i+l]);
    	mpz_clear(exponentation);
    }
    //Applies the w rounds over the two branches
    for (int i = 0; i <= w - 1 ; i++ ){
    	mpz_t exponentation1, exponentation2,sum2,sum1, multiplication1;
    	mpz_inits(exponentation1, exponentation2,sum2,sum1,multiplication1, NULL);

        if ( i % 2 == 0 ) {

        	mpz_ui_pow_ui (exponentation1, (unsigned long int) 2 , (unsigned long int) (f-32) );
    		mpz_ui_pow_ui (exponentation2, (unsigned long int) s , (unsigned long int) l );
    		mpz_mul_ui (multiplication1, exponentation1, (unsigned long int) (Tr ^ i) );
    		mpz_add ( sum1, multiplication1 , R[i] );
    		mpz_add ( sum2, L[i] , FK(sum1) );
    		mpz_mod (L[i+1], sum2, exponentation2);
    		mpz_set (R[i+1], R[i] );
            /*
            These lines are leaved commented for more clarity while the code is debugged
            L[i+1] = ( L[i] + FK((Tr ^ i)* square_and_multiply( 2 , f-32 ) + R[i]) ) % square_and_multiply( s , l );
            R[i+1] = R[i];*/
        }
        else{
        	mpz_ui_pow_ui (exponentation1, (unsigned long int) 2 , (unsigned long int) (f-32) );
    		mpz_ui_pow_ui (exponentation2, (unsigned long int) s , (unsigned long int) r );
    		mpz_mul_ui (multiplication1, exponentation1, (unsigned long int) (Tl ^ i) );
    		mpz_add ( sum1, multiplication1 , L[i] );
    		mpz_add ( sum2, R[i] , FK(sum1) );
    		mpz_mod (R[i+1], sum2, exponentation2);
    		mpz_set (L[i+1], L[i] );
    		/* These lines are leaved commented for more clarity while the code is debugged
            R[i+1] = ( R[i] + FK((Tl ^ i)* square_and_multiply( 2 , f-32 ) + L[i]) ) % square_and_multiply( s , r );
            L[i+1] = L[i]; */
        }
        mpz_clears(exponentation1, exponentation2,sum2,sum1,multiplication1, NULL);

    }
    //Transforms back to the s-base 
    for (int i = 0 ; i <= l - 1 ; i ++  ){
    	mpz_t residue,intermediate;
    	mpz_inits(residue,intermediate, NULL);
    	mpz_mod_ui (residue, L[w], (unsigned long int) s);
    	Y[i] = mpz_get_si (residue);
    	
    	//TODO: I'm thinking is enough an exact division
    	mpz_sub (intermediate, L[w], residue);
    	mpz_cdiv_q_ui  (L[w], intermediate, (unsigned long int) s);
    	/* These lines are leaved commented for more clarity while the code is debugged
        Y[i] = L[w] % s;
        L[w] = ( L[w] - Y[i] ) / s;
        */
        mpz_clears(residue,intermediate, NULL);
    }
    for (int i = 0 ; i <= r - 1 ; i ++  ){
    	mpz_t residue,intermediate;
    	mpz_inits(residue,intermediate, NULL);
    	mpz_mod_ui (residue, R[w], (unsigned long int) s);
    	Y[i + l] = mpz_get_si (residue);
    	
    	//TODO: I'm thinking is enough an exact division
    	mpz_sub (intermediate, R[w], residue);
    	mpz_cdiv_q_ui  (R[w], intermediate, (unsigned long int) s);
    	/* These lines are leaved commented for more clarity while the code is debugged
        Y[ i + l ] = R[w] % s;
        R[w] = ( R[w] - Y[ i + l ] ) / s;
        */
        mpz_clears(residue,intermediate, NULL);
        
    }
    for (int i = 0 ; i < w + 1 ; i++ )
    	mpz_clears( L[i] , R[i] , NULL );
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
vector<long long> InternalBlockDecipher(/*TODO: research how to insert a function here*/long long s, long long b
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

static void encrypt_aes_()
{
    // Example of more verbose verification

    uint8_t i, buf[64], buf2[64];

    // 128bit key
    uint8_t key[16] =        { (uint8_t) 0x2b, (uint8_t) 0x7e, (uint8_t) 0x15, (uint8_t) 0x16, (uint8_t) 0x28, (uint8_t) 0xae, (uint8_t) 0xd2, (uint8_t) 0xa6, (uint8_t) 0xab, (uint8_t) 0xf7, (uint8_t) 0x15, (uint8_t) 0x88, (uint8_t) 0x09, (uint8_t) 0xcf, (uint8_t) 0x4f, (uint8_t) 0x3c };
    // 512bit text
    uint8_t plain_text[64] = { (uint8_t) 0x6b, (uint8_t) 0xc1, (uint8_t) 0xbe, (uint8_t) 0xe2, (uint8_t) 0x2e, (uint8_t) 0x40, (uint8_t) 0x9f, (uint8_t) 0x96, (uint8_t) 0xe9, (uint8_t) 0x3d, (uint8_t) 0x7e, (uint8_t) 0x11, (uint8_t) 0x73, (uint8_t) 0x93, (uint8_t) 0x17, (uint8_t) 0x2a,
                               (uint8_t) 0xae, (uint8_t) 0x2d, (uint8_t) 0x8a, (uint8_t) 0x57, (uint8_t) 0x1e, (uint8_t) 0x03, (uint8_t) 0xac, (uint8_t) 0x9c, (uint8_t) 0x9e, (uint8_t) 0xb7, (uint8_t) 0x6f, (uint8_t) 0xac, (uint8_t) 0x45, (uint8_t) 0xaf, (uint8_t) 0x8e, (uint8_t) 0x51,
                               (uint8_t) 0x30, (uint8_t) 0xc8, (uint8_t) 0x1c, (uint8_t) 0x46, (uint8_t) 0xa3, (uint8_t) 0x5c, (uint8_t) 0xe4, (uint8_t) 0x11, (uint8_t) 0xe5, (uint8_t) 0xfb, (uint8_t) 0xc1, (uint8_t) 0x19, (uint8_t) 0x1a, (uint8_t) 0x0a, (uint8_t) 0x52, (uint8_t) 0xef,
                               (uint8_t) 0xf6, (uint8_t) 0x9f, (uint8_t) 0x24, (uint8_t) 0x45, (uint8_t) 0xdf, (uint8_t) 0x4f, (uint8_t) 0x9b, (uint8_t) 0x17, (uint8_t) 0xad, (uint8_t) 0x2b, (uint8_t) 0x41, (uint8_t) 0x7b, (uint8_t) 0xe6, (uint8_t) 0x6c, (uint8_t) 0x37, (uint8_t) 0x10 };

    memset(buf, 0, 64);
    memset(buf2, 0, 64);

    // print text to encrypt, key and IV
    printf("ECB encrypt verbose:\n\n");
    printf("plain text:\n");
    for(i = (uint8_t) 0; i < (uint8_t) 4; ++i)
    {
        phex(plain_text + i * (uint8_t) 16);
    }
    printf("\n");

    printf("key:\n");
    phex(key);
    printf("\n");

    // print the resulting cipher as 4 x 16 byte strings
    printf("ciphertext:\n");
    for(i = 0; i < 4; ++i)
    {
        AES128_ECB_encrypt(plain_text + (i*16), key, buf+(i*16));
        phex(buf + (i*16));
    }
    printf("\n");
}

int main(){

    return 0;
}
