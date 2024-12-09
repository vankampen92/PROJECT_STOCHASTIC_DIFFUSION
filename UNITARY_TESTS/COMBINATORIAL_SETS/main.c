#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

int CombiNumber(int , int );
double logCombiNumber(double , int ); 
void GSL_Init_Random_Seed(const gsl_rng * );
void GSL_Init_Random_Seed_from_File(const gsl_rng * );
void Random_Combinatorial_Set (int n, int k, int * IX);

gsl_rng * r; /* Global generator defined in main.c */

// #define VERBOSE

#define CURRENT_TIME_RDN_SEED

int
main (void)
{
  int j, m_0, m, n, k, l; 
  gsl_combination * c;
  size_t i;
  size_t m_i; 
  int NEXT; 

  /* GNU Scientific Library:  Random numbers set up */
  const gsl_rng_type * T;
  const gsl_rng_type **t, **t0;

  gsl_rng_env_setup();  /* Environment variables  
                          GSL_RNG_TYPE and GSL_RNG_SEED are read  
                          to set the corresponding library variables 
                          gsl_rng_default and gsl_rng_default_seed. 
                        */
  T = gsl_rng_taus2;    //T = gsl_rng_default; T = gsl_rng_ranlux389;
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

  /* End of Random number setting */

  #if defined VERBOSE
  /* BEGIN: Checking Random Number Generator Setup */
  for(j=0; j<10; j++){
    printf( "f(%d)=%g, ", j, gsl_rng_uniform(r) );
    printf( "f_GAUS(%d)=%g\n", j, gsl_ran_gaussian(r, 1.0) );
  }
  printf("\n"); getchar();
  /*   END: Checking Random Number Generator Setup */
  #endif

  /* Combinations from n elements taken in subsets of size k */
  n = 5;
  k = 2;  
  m = CombiNumber(n, k); /* CombiNumber(...) (stat.c) */
  /*  unsigned long int gsl_rng_uniform_int(const gsl_rng *r, unsigned long int m)

      This function returns a random integer from 0 to m-1 inclusive by scaling down 
      and/or discarding samples from the generator r. All integers in the range [0,m-1] 
      are produced with equal probability.
  */
  printf(" Total number of combinations from n = %d elements taken in subsets of size k = %d ", 
           n, k);
  printf("is the combinatorial number (n over k): %d\n\n", m);

  m_0 = gsl_rng_uniform_int(r, m);      
  /* The m_0-th combination will be selected from the orderd sets of potential 
     combinations with equal probabiltiy, and then retrived in a vector. 
  */
  printf(" Selection of the %d-th combination\n", m_0);
  printf (" All subsets of {0, ..., %d} of %d elements:\n", n-1, k);
  i = (size_t)k; 

  c = gsl_combination_calloc (n, i);
  do
    {
      printf ("{");   
      gsl_combination_fprintf (stdout, c, " %u");
      printf (" }\n");
    }
  while (gsl_combination_next (c) == GSL_SUCCESS);
  gsl_combination_free (c);

  printf("\n");
  printf (" (1) All subsets of {0, ..., %d} of %d elements up to the %d-th:\n", n-1, k, m_0);
  c = gsl_combination_calloc (n, i);
  NEXT = 1; j = 0;   
  while (NEXT == 1 && j<= m_0) {

    printf ("{");
    gsl_combination_fprintf (stdout, c, " %u");
    printf (" }\n");

    if (gsl_combination_next (c) == GSL_SUCCESS) {
      NEXT = 1; 
      j++;
    }
    else 
      NEXT = 0;     
  }
  gsl_combination_free (c);

  printf("\n");
  printf (" (2) All subsets of {0, ..., %d} of %d elements up to the %d-th:\n", n-1, k, m_0);
  c = gsl_combination_calloc (n, i);
  int * IX = (int *)calloc(k, sizeof(int));
  j = 0;   
  do
    {
      printf ("{");
      gsl_combination_fprintf (stdout, c, " %u");
      printf (" }\n");

      for(l = 0; l<k; l++) {
        m_i = (size_t)l; 
        IX[l] = (int)gsl_combination_get(c, m_i);
      }

      j++; 
    }
  while (gsl_combination_next (c) == GSL_SUCCESS && j <= m_0);
  gsl_combination_free (c);

  printf("\n");
  printf(" The required %d-th combination is: \n", m_0);
  printf ("{");
    for(l = 0; l<k; l++) 
        printf(" %d", IX[l]);                        
  printf (" }\n");

  getchar(); 
  Random_Combinatorial_Set(n, k, IX);

  free(IX);

  return 0;
}

void Random_Combinatorial_Set (int n, int k, int * IX)
{
  /* This function returns the IX array which strores a random 
     combinatorial set amb all possible combinations of n elements, 
     when taken in subsets of size k. 
  */
  int l, j, m, m_0; 
  size_t i, m_i; 
  gsl_combination * c;

  m = CombiNumber(n, k); /* CombiNumber(...) (stat.c) */
  /*  unsigned long int gsl_rng_uniform_int(const gsl_rng *r, unsigned long int m)

      This function returns a random integer from 0 to m-1 inclusive by scaling down 
      and/or discarding samples from the generator r. All integers in the range [0,m-1] 
      are produced with equal probability.
  */
#if defined VERBOSE
  printf(" Total number of combinations from n = %d elements taken in subsets of size k = %d ", 
           n, k);
  printf("is the combinatorial number (n over k): %d\n\n", m);
#endif 

  m_0 = gsl_rng_uniform_int(r, m);      
  /* The m_0-th combination will be selected from the orderd sets of potential 
     combinations with equal probabiltiy, and then retrived in the vector IX[]. 
  */
  i = (size_t)k; 
  c = gsl_combination_calloc (n, i);
#if defined VERBOSE  
  printf("\n");
  printf (" (2) All subsets of {0, ..., %d} of %d elements (up to the %d-th):\n", n-1, k, m_0);
#endif

  j = 0;   
  do
    {
#if defined VERBOSE
      printf ("{");
      gsl_combination_fprintf (stdout, c, " %u");
      printf (" }\n");
#endif

      for(l = 0; l<k; l++) {
        m_i = (size_t)l; 
        IX[l] = (int)gsl_combination_get(c, m_i);
      }

      j++; 
    }
  while (gsl_combination_next (c) == GSL_SUCCESS && j <= m_0);
  gsl_combination_free (c);

  printf("\n\n");
  printf(" The randomly chosen %d-th final combination is: \n", m_0);
  printf ("{");
    for(l = 0; l<k; l++) 
        printf(" %d", IX[l]);                        
  printf (" }\n");
}  
  
int CombiNumber(int n, int k)
{
  int a;
  double logof_a, x; 

  if(n < k){
    printf("No Combinatorial Number Exists...\n");
    printf("n = %d\tk = %d\n",n, k);
    printf("n is less than k!!!!\n");
    exit(0);
  }
  else if(k == n)
    a = 1;
  else{
    x       = (double)n;
    logof_a = logCombiNumber(x, k);     
  }

  a = (int)exp(logof_a);

  return(a);
}

double logCombiNumber(double x, int y)
{
  double a;
  int k;

  if(x < y){
    printf("No Combinatorial Number Exists...\n");
    printf("x = %g\ty = %d\n",x,y);
    printf("x is less than y!!!!\n");
    exit(0);
  }
  else if((int)x == y)
    a = 0.;
  else{
    a = 0.;
    for(k=0; k<y; k++)
      a += log((x-(double)k)/(y-(double)k));
  }

  return(a);
}

void GSL_Init_Random_Seed(const gsl_rng * r)
{
  /* This function seed the GSL Random Generator r
     with a seed which is different for each initialization
     according to current computer time 
  */
        unsigned long int     seed;
        time_t  nowtime;
        struct  tm *preztime;

        time(&nowtime);
        preztime = localtime(&nowtime);
        seed = (int)((preztime->tm_sec+1)*(preztime->tm_min+1)*
                (preztime->tm_hour+1)*(preztime->tm_year)*(preztime->tm_year));
        if(seed%2==0) seed++;

        printf(" Random Number Seed: %lu\n", seed);

        gsl_rng_set(r, seed);
}

void GSL_Init_Random_Seed_from_File(const gsl_rng * r)
{
        unsigned long int     seed;

        //seed = 100;       
        /* A script have changed the environmetal variable
           GSL_RNG_SEED before the execution of the program
           starts. This value has been set to gsl_rng_default_seed
           in gsl_random_number_Setup.c. This setup is 
           always done at the start of code execution.
        */

        printf ("GSL_RNG_SEED = %lu\n", gsl_rng_default_seed);

        gsl_rng_set(r, gsl_rng_default_seed);

        printf ("first value = %lu\n", gsl_rng_get (r));
}
