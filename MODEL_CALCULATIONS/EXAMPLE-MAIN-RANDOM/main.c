#ifndef THE_BASIC_C_HEADER_FILES
#define THE_BASIC_C_HEADER_FILES
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#endif

#ifndef GSL_HEADER_FILES
#define GSL_HEADER_FILES
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_heapsort.h>
#endif

#include <GSL_stat.h>

#define CURRENT_TIME_RDN_SEED
#define VERBOSE

gsl_rng * r;  /* Global Variable r to communicate with random gerenerator funcions */

int main()
{
  int i, n;

  // #include <gsl_random_number_Setup.c>
  
  /* B E G I N : Random init */
  const gsl_rng_type * T;
  const gsl_rng_type **t, **t0;

  gsl_rng_env_setup();
  
  T = gsl_rng_taus2; //T = gsl_rng_default; T = gsl_rng_ranlux389;

  r = gsl_rng_alloc(T);

  t0 = gsl_rng_types_setup ();

#if defined VERBOSE
  printf("Random number information... Available generators:\n");
  for(t=t0; *t != 0; t++){
    printf("%s, ", (*t)->name);
  }
#endif

  printf("\n In this example, the random generator at work... is the '%s' generator\n\n",
         gsl_rng_name(r));

#if defined CURRENT_TIME_RDN_SEED
  GSL_Init_Random_Seed(r);
#else
  GSL_Init_Random_Seed_from_File(r);
#endif

  /* BEGIN: Checking Random Number Generator Setup */
  for(i=0; i<10; i++){
    printf( "f(%d)=%g, ", i, gsl_rng_uniform(r) );
    printf( "f_GAUS(%d)=%g\n", i, gsl_ran_gaussian(r, 1.0) );
  }
  printf("\n"); 
  /*   END: Checking Random Number Generator Setup */
  /* E N D : Init Random */
  
  double * probability = (double *)calloc(10, sizeof(double) );

  probability[0] = 0.1;
  probability[1] = 0.4;
  probability[2] = 0.6;
  probability[3] = 0.4;
  probability[4] = 0.9;
  probability[5] = 0.23;
  probability[6] = 0.1;
  probability[7] = 0.8;
  probability[8] = 0.0;
  probability[9] = 0.9;

  /* BEGIN: Checking Sampling Discrete Probability */
  for(i=0; i<20; i++) {
    
    n = Discrete_Sampling(probability, 10);
    
    printf("Event(%d): %d\n", i, n);
    
  }
  /*   END: Checking Sampling Discrete Probability */
  
  free(probability);

  printf("\nEnd of progam\n");
  return(0);
}
