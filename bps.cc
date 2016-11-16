#include <iostream>
#include <vector>
#include <cmath>
#define TWEAK_BASE 2LL
#define TWEAK_LEFT_LENGTH 32LL

using namespace std;

struct bignum{
	int len;
	int64 data[bignumlen];
	int64 &operator [](int x){ return(data[x]);}
	const int64 &operator [](int x)const { return(data[x]);}
	bignum (){
		memset(data,0,sizeof(data));
		len=0;
	}
	void clear(){
		for(int i=len;i>=1;--i)data[i]=0;
		len=0;
	}
	int check (const bignum &a,const bignum &b){
		if(a.len>b.len)return(0);
		if(b.len>a.len)return(1);
		for(int i=a.len;i>=1;--i){
			if(a.data[i]<b.data[i])return(1);
			if(b.data[i]<a.data[i])return(0);
		}
		return 2;
	}
	bool operator < (const bignum &b){ return(check(*this,b)==1);}
	bool operator > (const bignum &b){ return(check(*this,b)==0);}
	bool operator <=(const bignum &b){ return(check(*this,b)>=1);}
	bool operator >=(const bignum &b){ return(check(*this,b)%2==0);}
	bool operator !=(const bignum &b){ return(check(*this,b)!=2);}
	bool operator ==(const bignum &b){ return(check(*this,b)==2);}

	bignum operator=(const bignum &x){
		for(int i=x.len+1;i<=len;++i)data[i]=0;
		for(int i=1;i<=x.len;++i)data[i]=x.data[i];
		len=x.len;
		return *this;
	}
	bignum operator=(int64 x){
		for(int i=len;i>=0;--i)data[i]=0;
		len=0;
		while(x){
			data[++len]=x%base;
			x/=base;
		}
		return *this;
	}
	bignum(int64 x){
		memset(data,0,sizeof(data));
		len=0;
		(*this)=x;
	}
	bignum operator *(const bignum &b){
		int i,j;
		bignum tmp;
		for(i=1;i<=len;++i)if(data[i]!=0)
			for(j=1;j<=b.len;++j)if(b.data[j]!=0){
				tmp.data[i+j-1]+=data[i]*b.data[j];
				tmp.data[i+j]+=tmp.data[i+j-1]/base;
				tmp.data[i+j-1]%=base;
			}
		tmp.len=len+b.len-1;
		while(tmp.data[tmp.len+1])tmp.len++;
		return tmp;
	}
	bignum operator *(int64 x){
		int i;
		bignum tmp;
		for(i=1;i<=len;++i)tmp[i]=data[i]*x;
		tmp.len=len;
		for(i=1;i<=len;++i){
			tmp[i+1]+=tmp[i]/base,tmp[i]%=base;
			if(tmp[i+1]&&i+1>tmp.len)tmp.len++;
		}
		return tmp;
	}
	bignum operator /(int64 x){
		int i;
		bignum tmp;
		int64 y=0;
		for(i=len;i>=1;--i){
			y=y*base+data[i];
			tmp[i]=y/x;
			y%=x;
		}
		tmp.len=len;
		while(tmp[tmp.len]==0&&tmp.len>1)tmp.len--;
		return tmp;
	}
	bignum operator /(const bignum &b){
		if(b.len<=1 && b[1]==0){
			printf("error! 被0除!");
			for(;;);
		}
		int i,l1=(len-1)*Blen,l2=(b.len-1)*Blen;
		int64 x=data[len],y=b[b.len];
		while(x)x/=10,l1++;
		while(y)y/=10,l2++;
		bignum tmp,chu,B;
		chu=*this; B=b;

		for(i=1;i*Blen<=l1-l2;++i)B*=base;
		for(i=1;i<=(l1-l2)%Blen;++i)B*=10;
		for(i=l1-l2;i>=0;--i){
			x=0;
			while(chu>=B)chu-=B,x++;
			tmp[i/Blen+1]=tmp[i/Blen+1]*10+x;
			B/=10;
		}
		tmp.len=(l1-l2)/Blen+1;
		while(tmp.len>=1 && !tmp[tmp.len])tmp.len--;
		return tmp;
	}
	bignum operator +(const bignum &b){
		bignum tmp;
		int i,l=max(len,b.len);
		for(i=1;i<=l;++i)tmp[i]=data[i]+b[i];
		for(i=1;i<=l;++i)tmp[i+1]+=tmp[i]/base,tmp[i]%=base;
		tmp.len=l;
		if(tmp[tmp.len+1])tmp.len++;
		return tmp;
	}
	bignum operator +(int64 x){
		bignum tmp; tmp=*this;
		tmp[1]+=x;
		for(int i=1;i<=len&&tmp[i]>=base;++i)tmp[i+1]+=tmp[i]/base,tmp[i]%=base;
		while(tmp[tmp.len+1])tmp.len++;
		return tmp;
	}
	bignum operator -(const bignum &b){
		int i;
		bignum tmp;
		for(i=1;i<=len;++i)tmp.data[i]=data[i]-b.data[i];
		for(i=1;i<=len;++i){
			if(tmp[i]<0)tmp.data[i]+=base,tmp.data[i+1]--;
		}
		tmp.len=len;
		while(tmp[tmp.len]==0&&tmp.len>1)tmp.len--;
		return tmp;
	}
	bignum operator -(int64 x){
		bignum tmp; tmp=*this;
		tmp[1]-=x;
		for(int i=1;i<=len&&tmp[i]<0;++i){
			tmp[i+1]+=(tmp[i]+1)/base-1;
			tmp[i]=(tmp[i]+1)%base+base-1;
		}
		while(!tmp[tmp.len]&&tmp.len>1)tmp.len--;
		return tmp;
	}
	int64 operator %(int64 x){
		int i;
		int64 y=0;
		for(i=len;i>=1;--i)y=(y*base+data[i])%x;
		return y;
	}
	bignum operator %(const bignum &b){
		if(b.len<=1 && b[1]==0){
			printf("error! 被0 mod!");
			for(;;);
		}
		int i,l1=(len-1)*Blen,l2=(b.len-1)*Blen;
		int64 x=data[len],y=b[b.len];
		while(x)x/=10,l1++;
		while(y)y/=10,l2++;
		bignum chu,B;
		chu=*this; B=b;

		for(i=1;i*Blen<=l1-l2;++i)B*=base;
		for(i=1;i<=(l1-l2)%Blen;++i)B*=10;
		for(i=l1-l2;i>=0;--i){
			while(chu>=B)chu-=B;
			B/=10;
		}
		return chu;
	}

	bignum operator +=(const bignum &b){return *this=(*this+b);}
	bignum operator *=(const bignum &b){return *this=(*this*b);}
	bignum operator -=(const bignum &b){return *this=(*this -b);}
	bignum operator /=(const bignum &b){return *this=(*this/b);}
	bignum operator %=(const bignum &b){return *this=(*this%b);}
	bignum operator *=(int64 x) {return( *this=(*this *x));}
	bignum operator +=(int64 x) {return( *this=(*this +x));}
	bignum operator -=(int64 x) {return( *this=(*this -x));}
	bignum operator /=(int64 x) {return( *this=(*this /x));}
	void read(){
		char c[bignumlen*Blen+10];
		scanf("%s",c+1);
		int l=strlen(c+1);
		(*this).clear();
		int64 x;
		for(int i=1;i<=(l-1)/Blen+1;++i){
			x=0;
			for(int j=l-Blen*i+1;j<=l-Blen*i+Blen;++j)if(j>=1)x=x*10+c[j]-48;
			data[++len]=x;
		}
	}
	void write(){
		printf("%I64d",data[len]);
		for(int i=len-1;i>=1;--i)printf("%0*I64d",Blen,data[i]);
	}
};


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

    vector< bignum > L( w + 1 ,0) , R( w + 1 ,0), Y(b) ;
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
