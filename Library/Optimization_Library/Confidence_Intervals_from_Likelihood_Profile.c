#include <MODEL.h>
/* 
   Confidence intervals are calculated by the profile likelihood method. 
   1.92 corresponds to (1-alpha)-percentile of the Chi^2(1) distribution
   (with one degree of freedom) where alpha is a significance level of 
   0.05

   Usually, 1.92 is approximated by 2.0 and given as the LIKELIHOOD_JUMP
   input paremter of the function. 
*/

void Confidence_Intervals_from_Likelihood_Profile( double * Profile, double * X, int N,
						   double LIKELIHOOD_JUMP, 
						   double * X_MLE, double * CI )
{
  /* Input arguments:
     . Profile[]: Likelihood vaalues for each parameter value  
     . X[]:       Vector of parameter values for which a likelihood has been evaluated
     . N : Number of discretizing points fof the likelhood profile

     Output arguments:
      . ( CI[0], CI[1] ): Confidence Interaval, 
      . * X_MLE   m.l.e of the model parameter
  */ 
  int i, i_MLE, n_BOOL; 
  gsl_vector * PROFILE = gsl_vector_alloc( N );
  double P_MLE; 
  
  for (i=0; i<N; i++) gsl_vector_set(PROFILE, i, Profile[i]);

  i_MLE = gsl_vector_min_index ( PROFILE );

  if (i_MLE > 0 && i_MLE < N-1 ) {  // There is a true interior minimum  
    * X_MLE      = X[i_MLE];
    P_MLE = Profile[i_MLE];
  
    i=0; n_BOOL=0; 
    while ( i < (1+i_MLE) && n_BOOL == 0 ) {
      if( (Profile[i_MLE - i] - P_MLE) > LIKELIHOOD_JUMP ) { 
	n_BOOL = 1;
	CI[0] = X[i_MLE - i];
      }
      i++;
      
      if( (i_MLE-i) == -1 ) break; 
    }
    
    i=0; n_BOOL=0; 
    while ( i < (N-i_MLE) && n_BOOL == 0 ) {
      if( (Profile[i_MLE + i] - P_MLE) > LIKELIHOOD_JUMP ) { 
	n_BOOL = 1;
	CI[1] = X[i_MLE + i];
      }
      i++;
      
      if( (i_MLE+i) == N ) break; 
    }
  }
  
  printf(" --------------------------------------------------------------------\n");
  printf( " m.l.e. parameter value from likelihood profile = %g\n", * X_MLE);
  printf( " Confidence Interaval:\n");
  printf( " [ %g,  (value = %g),  %g ]\n", CI[0], * X_MLE, CI[1] );
  printf(" --------------------------------------------------------------------\n");

  gsl_vector_free(PROFILE);
}
