#include <iostream>
#include <stdarg.h>
#include <obstack.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <stdint.h>
#include <algorithm>
#include "tiny-AES128-C/aes.c"

#define TWEAK_BASE 2LL
#define TWEAK_LEFT_LENGTH 32LL
#define F 128 //AES_BLOCK_LENGTH

using namespace std;


// prints string as hex
static void phex(uint8_t* str)
{
    unsigned char i;
    for(i = 0; i < 16; ++i)
        printf("%.2x", str[i]);
    printf("\n");
}


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


void encrypt_aes( mpz_t ciphertext  , const mpz_t message, const int K[] )
{
    // Example of more verbose verification

    uint8_t buf[16];
    uint8_t key[16];
    uint8_t plain_text[64];
    mpz_t q, r, temp ;
    mpz_inits( q , r , temp , NULL);

    for(int i = 0; i < 16; i++){
		key[i] = (uint8_t) K[i];
	}
    
    
    mpz_set( temp, message );
    for (int i = 15 ; i >= 0 ; i-- ){
    	mpz_fdiv_qr_ui ( q , r , temp, (unsigned long int) 256) ;
    	mpz_set (temp, q );
    	plain_text[i] = (uint8_t) mpz_get_ui (r);
    }

    /*

    // 128bit key
    uint8_t key[16] =        { (uint8_t) 0x2b, (uint8_t) 0x7e, (uint8_t) 0x15, (uint8_t) 0x16, (uint8_t) 0x28, (uint8_t) 0xae, (uint8_t) 0xd2, (uint8_t) 0xa6, (uint8_t) 0xab, (uint8_t) 0xf7, (uint8_t) 0x15, (uint8_t) 0x88, (uint8_t) 0x09, (uint8_t) 0xcf, (uint8_t) 0x4f, (uint8_t) 0x3c };
    // 512bit text
    uint8_t plain_text[64] = { (uint8_t) 0x6b, (uint8_t) 0xc1, (uint8_t) 0xbe, (uint8_t) 0xe2, (uint8_t) 0x2e, (uint8_t) 0x40, (uint8_t) 0x9f, (uint8_t) 0x96, (uint8_t) 0xe9, (uint8_t) 0x3d, (uint8_t) 0x7e, (uint8_t) 0x11, (uint8_t) 0x73, (uint8_t) 0x93, (uint8_t) 0x17, (uint8_t) 0x2a,
                               (uint8_t) 0xae, (uint8_t) 0x2d, (uint8_t) 0x8a, (uint8_t) 0x57, (uint8_t) 0x1e, (uint8_t) 0x03, (uint8_t) 0xac, (uint8_t) 0x9c, (uint8_t) 0x9e, (uint8_t) 0xb7, (uint8_t) 0x6f, (uint8_t) 0xac, (uint8_t) 0x45, (uint8_t) 0xaf, (uint8_t) 0x8e, (uint8_t) 0x51,
                               (uint8_t) 0x30, (uint8_t) 0xc8, (uint8_t) 0x1c, (uint8_t) 0x46, (uint8_t) 0xa3, (uint8_t) 0x5c, (uint8_t) 0xe4, (uint8_t) 0x11, (uint8_t) 0xe5, (uint8_t) 0xfb, (uint8_t) 0xc1, (uint8_t) 0x19, (uint8_t) 0x1a, (uint8_t) 0x0a, (uint8_t) 0x52, (uint8_t) 0xef,
                               (uint8_t) 0xf6, (uint8_t) 0x9f, (uint8_t) 0x24, (uint8_t) 0x45, (uint8_t) 0xdf, (uint8_t) 0x4f, (uint8_t) 0x9b, (uint8_t) 0x17, (uint8_t) 0xad, (uint8_t) 0x2b, (uint8_t) 0x41, (uint8_t) 0x7b, (uint8_t) 0xe6, (uint8_t) 0x6c, (uint8_t) 0x37, (uint8_t) 0x10 };
    */

    memset(buf, 0, 16);


    AES128_ECB_encrypt( plain_text , key, buf );

    for (int i = 0; i < 16; ++i){
    	mpz_mul_ui (temp, ciphertext , 256);
    	mpz_add_ui ( ciphertext , temp , (unsigned long int) buf[i] );
    }
    mpz_clears( q , r , temp , NULL );
    
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
vector<int> InternalBlockCipher(/*TODO: research how to insert a function here*/int s, int b
                        , int w, int X[], int K[], unsigned long long T){

    long long Tr = T % square_and_multiply( TWEAK_BASE , TWEAK_LEFT_LENGTH );

    long long Tl = (T - Tr) / square_and_multiply( TWEAK_BASE , TWEAK_LEFT_LENGTH);


      


    long long l = ceil(b/2.0), r = floor( b/2.0 );



    //Declare the BigIntegers to be used in each round 
    mpz_t  *L = new mpz_t[ w + 1 ]; 
    mpz_t  *R = new mpz_t[ w + 1 ]; 
    

    //Initialize the BigIntegers nums
    for (int i = 0 ; i < w + 1 ; i++ )
    	mpz_inits( L[i] , R[i] , NULL );

    vector<int> Y(b);

    //Transforms left branch from s-base to decimal base
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
    	mpz_t exponentation1, exponentation2,sum2,sum1, multiplication1,ciphertext;
    	mpz_inits(exponentation1, exponentation2,sum2,sum1,multiplication1,ciphertext, NULL);
    
        if ( i % 2 == 0 ) {

        	mpz_ui_pow_ui (exponentation1, (unsigned long int) 2 , (unsigned long int) (F-32) );
    		mpz_ui_pow_ui (exponentation2, (unsigned long int) s , (unsigned long int) l );
    		mpz_mul_ui (multiplication1, exponentation1, (unsigned long int) (Tr ^ i) );
    		mpz_add ( sum1, multiplication1 , R[i] );
     
    		encrypt_aes( ciphertext, sum1 , K );
     
    		mpz_add ( sum2, L[i] , ciphertext );     
    		mpz_mod (L[i+1], sum2, exponentation2);
    		mpz_set (R[i+1], R[i] );
            /*
            These lines are leaved commented for more clarity while the code is debugged
            L[i+1] = ( L[i] + FK((Tr ^ i)* square_and_multiply( 2 , f-32 ) + R[i]) ) % square_and_multiply( s , l );
            R[i+1] = R[i];*/
        }
        else{
        	mpz_ui_pow_ui (exponentation1, (unsigned long int) 2 , (unsigned long int) (F-32) );
    		mpz_ui_pow_ui (exponentation2, (unsigned long int) s , (unsigned long int) r );
    		mpz_mul_ui (multiplication1, exponentation1, (unsigned long int) (Tl ^ i) );
    		mpz_add ( sum1, multiplication1 , L[i] );
  
    		encrypt_aes( ciphertext, sum1 , K );
  
    		mpz_add ( sum2, R[i] , ciphertext );
    		mpz_mod (R[i+1], sum2, exponentation2);
    		mpz_set (L[i+1], L[i] );
    		/* These lines are leaved commented for more clarity while the code is debugged
            R[i+1] = ( R[i] + FK((Tl ^ i)* square_and_multiply( 2 , f-32 ) + L[i]) ) % square_and_multiply( s , r );
            L[i+1] = L[i]; */
        }
        mpz_clears(exponentation1, exponentation2,sum2,sum1,multiplication1,ciphertext, NULL);
  
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
vector<int> InternalBlockDecipher(/*TODO: research how to insert a function here*/int s, int b
                        , int w, int X[], int K[], unsigned long long T){

    long long Tr = T % square_and_multiply( TWEAK_BASE , TWEAK_LEFT_LENGTH );
    long long Tl = (T - Tr) / square_and_multiply( TWEAK_BASE , TWEAK_LEFT_LENGTH);

    

    long long l = ceil(b/2.0), r = floor( b/2.0 );


    //Declare the BigIntegers to be used in each round 
    mpz_t  *L = new mpz_t[ w + 1 ]; 
    mpz_t  *R = new mpz_t[ w + 1 ]; 
    

    //Initialize the BigIntegers nums
    for (int i = 0 ; i < w + 1 ; i++ )
        mpz_inits( L[i] , R[i] , NULL );

    vector<int> Y(b);

    //Transforms left branch from s-base to decimal base
    for (int i = 0 ; i < l ; i++ ){
        mpz_t exponentation;
        mpz_init(exponentation); 
        mpz_ui_pow_ui (exponentation, (unsigned long int) s , (unsigned long int) i );
        mpz_addmul_ui (L[w], exponentation, (unsigned long int) X[i]);
        mpz_clear(exponentation);
    }



    //Transforms right branch from s-base to decimal base
    for (int i = 0 ; i < r ; i++ ){
        mpz_t exponentation;
        mpz_init(exponentation); 
        mpz_ui_pow_ui (exponentation, (unsigned long int) s , (unsigned long int) i );
        mpz_addmul_ui (R[w], exponentation, (unsigned long int) X[i+l]);
        mpz_clear(exponentation);
    }
    //Applies the w rounds over the two branches
    for (int i = w-1; i >= 0 ; i-- ){
        mpz_t exponentation1, exponentation2,sum2,sum1, multiplication1,ciphertext;
        mpz_inits(exponentation1, exponentation2,sum2,sum1,multiplication1,ciphertext, NULL);
     
        if ( i % 2 == 0 ) {
            mpz_ui_pow_ui (exponentation1, (unsigned long int) 2 , (unsigned long int) (F-32) );
            mpz_ui_pow_ui (exponentation2, (unsigned long int) s , (unsigned long int) l );
            mpz_mul_ui (multiplication1, exponentation1, (unsigned long int) (Tr ^ i) );
            mpz_add ( sum1, multiplication1 , R[i+1] );     

     
            encrypt_aes( ciphertext, sum1 , K );
     
            mpz_sub ( sum2, L[i+1] , ciphertext );
     
            mpz_mod (L[i], sum2, exponentation2);
            mpz_set (R[i], R[i+1] );
            /*
            These lines are leaved commented for more clarity while the code is debugged
            L[i+1] = ( L[i] + FK((Tr ^ i)* square_and_multiply( 2 , f-32 ) + R[i]) ) % square_and_multiply( s , l );
            R[i+1] = R[i];*/
        }
        else{
            mpz_ui_pow_ui (exponentation1, (unsigned long int) 2 , (unsigned long int) (F-32) );
            mpz_ui_pow_ui (exponentation2, (unsigned long int) s , (unsigned long int) r );
            mpz_mul_ui (multiplication1, exponentation1, (unsigned long int) (Tl ^ i) );
            mpz_add ( sum1, multiplication1 , L[i+1] );
     
          
       
            encrypt_aes( ciphertext, sum1 , K );
       
            mpz_sub ( sum2, R[i+1] , ciphertext );
            mpz_mod (R[i], sum2, exponentation2);
            mpz_set (L[i], L[i+1] );
            /* These lines are leaved commented for more clarity while the code is debugged
            R[i+1] = ( R[i] + FK((Tl ^ i)* square_and_multiply( 2 , f-32 ) + L[i]) ) % square_and_multiply( s , r );
            L[i+1] = L[i]; */
        }
        mpz_clears(exponentation1, exponentation2,sum2,sum1,multiplication1,ciphertext, NULL);     
    }
    //Transforms back to the s-base 
    for (int i = 0 ; i <= l - 1 ; i ++  ){
        mpz_t residue,intermediate;
        mpz_inits(residue,intermediate, NULL);
        mpz_mod_ui (residue, L[0], (unsigned long int) s);
        Y[i] = mpz_get_si (residue);
        
        //TODO: I'm thinking is enough an exact division
        mpz_sub (intermediate, L[0], residue);
        mpz_cdiv_q_ui  (L[0], intermediate, (unsigned long int) s);
        /* These lines are leaved commented for more clarity while the code is debugged
        Y[i] = L[w] % s;
        L[w] = ( L[w] - Y[i] ) / s;
        */
        mpz_clears(residue,intermediate, NULL);
    }
    for (int i = 0 ; i <= r - 1 ; i ++  ){
        mpz_t residue,intermediate;
        mpz_inits(residue,intermediate, NULL);
        mpz_mod_ui (residue, R[0], (unsigned long int) s);
        Y[i + l] = mpz_get_si (residue);
        
        //TODO: I'm thinking is enough an exact division
        mpz_sub (intermediate, R[0], residue);
        mpz_cdiv_q_ui  (R[0], intermediate, (unsigned long int) s);
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



int main(){
	int K[16] = {  0xEF,  0x43,  0x59,  0xD8,  0xD5,  0x80,  0xAA,  0x4F,  0x7F,  0x03,  0x6D,
					0x6F,  0x04,  0xFC,  0x6A,  0x94 };
	int X[18] = {8,9,0,1,2,1,2,3,4,5,6,7,8,9,0,0,0,0};
	int f = 128, s = 10 , b = 18 , w = 8 ;
	unsigned long long T = 0xD8E7920AFA330A73;

	vector<int> ciphertext =  InternalBlockCipher(s, b, w, X, K, T);
	cout << "The ciphertext is : "; 
	for (int i = 0; i < (int) ciphertext.size() ; ++i){
		cout << ciphertext[i] << " ";
	}
    vector<int> plaintext =  InternalBlockDecipher(s, b, w, &ciphertext[0], K, T);
    cout<< endl << "The plaintext is : "; 
    for (int i = 0; i < (int) plaintext.size() ; ++i){
        cout << plaintext[i] << " ";
    }
    cout<<endl;
    return 0;
}
