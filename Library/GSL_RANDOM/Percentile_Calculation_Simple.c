#include "GSL_stat.h"

void Percentile_Calculation_Simple(double *X, double ** Variable_0, 
				   int No, int No_REPLICATES, 
				   double * Per_Cent, int NoPer_Cent, 
				   double ** Y) 
{
  /* Input: 
   *
   *  . Variable_0 is an array of length No x No_REPLICATES:
   *  . For instance:
   *  .                Variable_0[i][j]
   *  .                           i = 0, ... , No_REPLICATES-1;
   *  .                           j = 0, ... , No-1 
   *  . Percentile vector, for instance, if NoPer_Cent = 5, 
   *                Per_Cent[5]  = { 0.05,  0.25,  0.5,  0.75,  0.95 } 
   *  . X is an array of length No 
   *
   * Output:
   *
   *  . Y[][] contains the percentile values in the following structure:
   *
   *                             X     p(0.05)   p(0.25)   p(0.5)   p(0.75)   p(0.95)
   *                             ...   ...       ...       ...      ...       ...
   *                             X[j]  Y[0][j]   Y[1][j]   Y[2][j]  Y[3][j]   Y[4][j]
   */
  int No_of_DATA;
  int i, j, k, n, no_per;
  double value;
  
  gsl_vector * xx;
  gsl_vector * x = gsl_vector_alloc(No_REPLICATES);
  gsl_permutation * p;  
  
  for (j = 0; j < No; j++){
    
    No_of_DATA = 0;
    for( n = 0; n < No_REPLICATES; n++ ) {
    
	    value = Variable_0[n][j];
	    gsl_vector_set(x, No_of_DATA, value);
	    
	    No_of_DATA++;
    }
        
    if (No_of_DATA > 0) {
      
      xx = gsl_vector_alloc(No_of_DATA);
      p = gsl_permutation_alloc(No_of_DATA);
      
      for( n=0; n<No_of_DATA; n++ ) {
	value = gsl_vector_get(x, n);
	gsl_vector_set(xx, n, value);
      }
      /* Here, array xx[] now is identical to x[] */
  
      gsl_sort_vector_index (p, xx);      
      
      for(k = 0; k < NoPer_Cent; k++){
	
	no_per  = Per_Cent[k] * No_of_DATA;

	assert(no_per > 0 && no_per < No_of_DATA);

	i       = gsl_permutation_get(p, no_per-1);
	
	value = gsl_vector_get( x, i );

	Y[k][j] = value;   
      }  

      gsl_permutation_free(p);
      gsl_vector_free(xx);
    }

    else {
      
      printf(" Number of Replicates is zero!\n");
      printf(" Percentile calculation is not possible!\n");
      printf(" No of REPLICATES: %d\n", No_REPLICATES);
      exit(0);
    }    
  }
  
  gsl_vector_free(x);
}
