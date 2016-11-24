
#include <stdarg.h>
#include <obstack.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <vector>


using namespace std;

int main(void) {
  mpz_t  testing[10] ;
  for (int i = 0 ; i < 10 ; i++ )
      mpz_init( testing[i]  );
  mpz_t result, base;
  mpz_inits(result,base,NULL);
  mpz_set_str(base, "2", 10);
  mpz_pow_ui(result, base, 64);
  mpz_out_str(stdout, 10, result);
  return 0;
}

