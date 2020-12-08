#include <MODEL.h>

double GSL_neglog_Error_Probability_Model( double * Data, double * Theory,
					   int N , int No_of_VARIABLES,
					   Parameter_Fitting  * F, 
	                                   double (* Error_Model )(double, double, Parameter_Fitting *) )
{
  int j, n;
  double * Value = (double *)calloc(N, sizeof(double) );
  double neg_X, Max_Value; 

  Parameter_Table * Table = F->Table;
  
  n = 0;
  for( j=1; j<N; j++ ) { 
    Value[j] = (* Error_Model)( Data[j], Theory[j], F );

    if (Value[j] <= 0.0) n++;

    if(F->Verbose == 1)
      printf("[n = %d]: NLL (t=%g\t Data: %g\t Theory: %g) = %g\n", n,
	     Table->T->Time_Vector[j], Data[j], Theory[j], Value[j]);
  }

  if(n > 1 && n < (N-1) ) {
    Max_Value = Value[1];
    for( j=1; j<N; j++ ) Max_Value = MAX(Max_Value, Value[j]);

    for( j=1; j<N; j++ )
      if( Value[j] <= 0.0 ) Value[j] = Max_Value;
  }

  if (n == (N-1)) {
    neg_X = DBL_MAX/(double)No_of_VARIABLES;
  }
  else {
    neg_X = 0.0;
    for( j=1; j<N; j++ ) neg_X += Value[j];
  }

  free(Value);
  
  return(neg_X); 
}

