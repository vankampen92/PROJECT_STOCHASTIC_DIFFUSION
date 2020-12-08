#include <MODEL.h>

int Checking_for_Parameter_Boundaries( Parameter_Fitting * F, const gsl_vector * x )
{
  /* Model parameters are correctly bounded if PAR_DEFINITION 
     flag is left unchanged and returned as 1. 
  */
  int i; 
  int PAR_DEFINITION;
  int No_of_PARAMETERS_0, No_of_PARAMETERS_1, No_of_PARAMETERS_2;
  int No_of_PARAMETERS; 
  double x_Value; 
  double x_Value_min; 
  double x_Value_MAX;
  
  No_of_PARAMETERS_0 = F->Space->No_of_PARAMETERS ;
  No_of_PARAMETERS_1 = F->Table->No_of_IC;
  No_of_PARAMETERS_2 = F->Table->No_of_ERROR_PARAMETERS;
  
  No_of_PARAMETERS = No_of_PARAMETERS_0 + No_of_PARAMETERS_1 + No_of_PARAMETERS_2; 
  
  PAR_DEFINITION = 1;  
  
  if( No_of_PARAMETERS != x->size ) {
    printf( "No_of_PARAMETERS = %d\t Subspace Dimension = %zu\n",
	    No_of_PARAMETERS, x->size );
    printf( "You are very much in troble. These two number should match!!!\n");
    printf( "The program will exit\n");
    exit(0);
  }

  for( i=0; i<No_of_PARAMETERS_0; i++ ){
    x_Value = gsl_vector_get( x, i);
    x_Value_min = gsl_vector_get(F->Space->P_min, i); 
    x_Value_MAX = gsl_vector_get(F->Space->P_MAX, i);
    
    if( x_Value < x_Value_min || x_Value > x_Value_MAX ) { 
      PAR_DEFINITION = 0;
      break;
    }
  }

  for( i=0; i<No_of_PARAMETERS_1; i++ ){
    x_Value = gsl_vector_get( x, i+No_of_PARAMETERS_0);
    x_Value_min = gsl_vector_get(F->Table->IC_Space->P_min, i); 
    x_Value_MAX = gsl_vector_get(F->Table->IC_Space->P_MAX, i);
    
    if( x_Value < x_Value_min || x_Value > x_Value_MAX ) { 
      PAR_DEFINITION = 0;
      break;
    }
  }

  for( i=0; i<No_of_PARAMETERS_2; i++ ){
    x_Value = gsl_vector_get( x, i+No_of_PARAMETERS_0+No_of_PARAMETERS_1);
    x_Value_min = gsl_vector_get(F->Table->E_Space->P_min, i); 
    x_Value_MAX = gsl_vector_get(F->Table->E_Space->P_MAX, i);
    
    if( x_Value < x_Value_min || x_Value > x_Value_MAX ) { 
      PAR_DEFINITION = 0;
      break;
    }
  }
  
  return (PAR_DEFINITION);
}
