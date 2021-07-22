#include <MODEL.h>

extern int MY_ERROR_HANDLER; 

void copy_accuracies_into_gsl_vector(gsl_vector * ss,
				     const gsl_vector * Space_Accuracy, int No_of_PARAMETERS,
				     const gsl_vector * IC_Space_Accuracy, int No_of_IC, 
				     const gsl_vector * E_Space_Accuracy, int No_of_ERROR_PARAMETERS)
{
  int i;
  int P_0, P_1, P_2;
  double value;
  
  P_0 = No_of_PARAMETERS;
  P_1 = No_of_IC;
  P_2 = No_of_ERROR_PARAMETERS;

  for( i = 0;  i < P_0; i++ ) {
    value = gsl_vector_get( Space_Accuracy, i);
    gsl_vector_set(ss, i, value);
  }
  for( i = P_0; i < (P_0+P_1) ; i++ ) {
    value = gsl_vector_get( IC_Space_Accuracy, i-P_0);
    gsl_vector_set(ss, i, value);
  }
  for( i = (P_0+P_1); i < (P_0+P_1+P_2); i++ ) {
    value = gsl_vector_get (E_Space_Accuracy, i-P_0-P_1);
    gsl_vector_set(ss, i, value);
  }
}

double GSL_Minimization_Simplex (Parameter_Fitting * F,
				 gsl_vector * Initial_Guess,
				 gsl_vector * Solution,
				 double ( * Function )( const gsl_vector * , void * ) )
{
  int i, key;
 
  double value;
  
  Parameter_Table * Table  = F->Table;
  Parameter_Space * Space  = F->Space;
  // Parameter_Model * P      = F->P;

  const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex;
  
  gsl_multimin_fminimizer * s = NULL;
  gsl_vector * ss;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  int P_0, P_1, P_2, No_of_PARAMETERS;
  
  P_0 = Space->No_of_PARAMETERS;
  P_1 = Table->No_of_IC;
  P_2 = Table->No_of_ERROR_PARAMETERS;
  
  No_of_PARAMETERS = P_0 + P_1 + P_2;

  double TOLERANCE         = Space->TOLERANCE;
  int MAX_No_of_ITERATIONS = Space->MAX_No_of_ITERATIONS;

  /* Set initial step sizes */
  ss = gsl_vector_alloc ( No_of_PARAMETERS );
  
  // gsl_vector_set_all (ss, 0.01);
  
  copy_accuracies_into_gsl_vector(ss,
				  Space->Accuracy, Space->No_of_PARAMETERS,
				  Table->IC_Space->Accuracy, Table->No_of_IC,
				  Table->E_Space->Accuracy, Table->No_of_ERROR_PARAMETERS);

  if ( F->Verbose == 1 ) printf("Total Parameter Space Dimension = %d\n", No_of_PARAMETERS);

  minex_func.n = No_of_PARAMETERS;
  minex_func.f = Function;
  minex_func.params = (void *)F;
  
  s = gsl_multimin_fminimizer_alloc (T, No_of_PARAMETERS);
  gsl_multimin_fminimizer_set (s, &minex_func, Initial_Guess, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      
      if (status) break;
      
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, TOLERANCE);

      if (status == GSL_SUCCESS)
	printf("Solution converged to a minimum!!!\n");

      if ( Checking_for_Parameter_Boundaries(F, s->x) == 0 ) {
	printf("The search will be aborted. Search is out of bounds!!!\n");
	break;
      }
      
      if (F->Verbose == 1) {
	
	if (status == GSL_SUCCESS) printf("converged to minimum at\n");
	
	printf("%s (%3d): ", F->Data->Name, (int)iter);
	for( i = 0;  i < P_0; i++ ) {
	  key = Space->Parameter_Index[i];
	  printf("%s=%10.3e   ", Table->Symbol_Parameters[key], gsl_vector_get (s->x, i) );
	}
	for( i = P_0; i < (P_0+P_1) ; i++ ) {
	  key = Table->IC_Space->Parameter_Index[i-P_0];
	  printf("%s=%10.3e   ", Table->Model_Variable_Symbol[key], gsl_vector_get (s->x, i) );
	}
	for( i = (P_0+P_1); i < (P_0+P_1+P_2); i++ ) {
	  key = Table->E_Space->Parameter_Index[i-P_0-P_1];
	  if(key < Table->OUTPUT_VARIABLES_GENUINE)
	    printf("Err(%s)=%10.3e   ",
		   Table->Output_Variable_Symbol[key], gsl_vector_get (s->x, i) );
	  else
	    printf("Error =%10.3e   ", gsl_vector_get (s->x, i) );
	}
	printf("\n");
	if (s->fval < DBL_MAX/(double)F->Data->No_of_POINTS) 
	  printf("NLL(Paramereters |Data) = %7.3f size = %.3f\n", s->fval, size);
	else 
	  printf("NLL(Paramereters |Data) = %s size = %.3f\n", "infty", size);

#if defined CPGPLOT_REPRESENTATION
	// GSL_CPGPLOT_Minimization_Simplex(F, s->x, iter, Function); 
	// getchar(); 
#endif
      } 
    }
  while (status == GSL_CONTINUE && iter < MAX_No_of_ITERATIONS );

  if (MY_ERROR_HANDLER == 0) {
#if defined CPGPLOT_REPRESENTATION
    GSL_CPGPLOT_Minimization_Simplex(F, s->x, iter, Function); 
#endif   
  }
  
  gsl_vector_memcpy (Solution, s->x );
  value = s->fval; // Min Value !!!
  
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  //return status;
  return( value );
}
