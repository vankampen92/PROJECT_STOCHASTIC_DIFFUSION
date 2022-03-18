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

gsl_rng * r;

void Update_Total_Rate(int K, double b, double dx, double a, double v,
		       double by, double dy, double * prob, double * R,
		       int n, int m, int l);

int main(){

  int i, i0;
  int s;
  
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

 
  double time, time_0;
  double b, dx, a, v, by, dy;
  int K;
  
  K = 2000; 
  b = 0.1;
  dx = 0.1;
  a = 0.1;
  v = 0.1;
  by = 0.1;
  dy = 0.1;
    
  int n, m, l, n0, m0, l0; 
  
  /* Initial Conditions */
  n0 = 1000;
  m0 = 1000;
  l0 = 0;
  time_0 = 0.0;
 
  double p0, p1, p2, p3, p4, p5, p6;
  double R;
  double * prob = (double *)calloc(7, sizeof(double) );
  
  p0 = b*(1 - n/(double)K)*(double)n;
  p1 = dx*(double)n;
  p2 = a*(double)n /K *(double)m;
  p3 = v*(double)l;
  p4 = by*(double)l;
  p5 = dy*(double)m;
  p6 = dy*(double)l;

  R = p0 + p1 + p2 + p3 + p4 + p5 + p6;

  prob[0] = p0/R;
  prob[1] = p1/R;
  prob[2] = p2/R;
  prob[3] = p3/R;
  prob[4] = p4/R;
  prob[5] = p5/R;
  prob[6] = p6/R;

  /* BEGIN: Checking Discrete Sampling */
  for(i=0; i<20; i++) {
    
    s = Discrete_Sampling(prob, 7);
    
    printf("Event(%d): %d\n", i, s);
  }
  /* E N D : */
  printf("Initiating Stochastic Temporal Evolution...\n"); 
  getchar(); /* Press any key to continue... */

  printf("Initial configuration at initial time:\n");

  n = n0;
  m = m0;
  l = l0;
  time = time_0;
  
  printf("Time = %g\t", time); 
  printf("Resource n = %d\t", n);
  printf("Serchers m = %d\t", m);
  printf("Handlers l = %d\n", l);
  printf("\n"); 

  for(i = 0; i < 100; i++){
  
    time += (-1/R * log ( gsl_rng_uniform_pos(r) )); 
      
    s = Discrete_Sampling(prob, 7) - 1;

    switch(s){
            case 0 : n = n + 1 ;
            break;
            case 1 : n = n - 1 ;
            break;
            case 2 : m = m - 1 ; l = l + 1;
            break;
            case 3 : m = m + 1; l = l - 1;
            break;
            case 4: m = m + 1;
            break;
            case 5: m = m - 1;
            break;
            case 6: l = l - 1;
            break;
    }
	
    printf("Time = %g\t", time); 
    printf("Resource n = %d\t", n);
    printf("Serchers m = %d\t", m);
    printf("Handlers l = %d\n", l);
    printf("\n"); 
    
    Update_Total_Rate(K, b, dx, a, v, by, dy, prob, &R,
		      n, m, l); 
  }
    
  free(prob);
    
  return 0;
}

void Update_Total_Rate(int K, double b, double dx, double a, double v,
		       double by, double dy, double * prob, double * R,
		       int n, int m, int l)
{		       
  double N, M, L;

  N = (double)n;
  M = (double)m;
  L = (double)l;
  
  double p0, p1, p2, p3, p4, p5, p6;
  
  p0 = b*(1-N/(double)K)*N;
  p1 = dx*N;
  p2 = a*(N/(double)K)*M;
  p3 = v*L;
  p4 = by*L;
  p5 = dy*M;
  p6 = dy*L;

  * R = p0 + p1 + p2 + p3 + p4 + p5 + p6;

  prob[0] = p0/(*R);
  prob[1] = p1/(*R);
  prob[2] = p2/(*R);
  prob[3] = p3/(*R);
  prob[4] = p4/(*R);
  prob[5] = p5/(*R);
  prob[6] = p6/(*R);
}
  
