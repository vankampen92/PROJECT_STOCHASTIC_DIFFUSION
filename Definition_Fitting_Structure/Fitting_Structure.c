#include <MODEL.h>

void Parameter_Fitting_Alloc( Parameter_Fitting * F, int Realizations, Parameter_Table * Table )  
{
  int i; 
  int No_of_MODEL_PARAMETERS = Table->TOTAL_No_of_MODEL_PARAMETERS;
  int No_of_IC               = Table->No_of_IC;
  int No_of_ERROR_PARAMETERS = Table->No_of_ERROR_PARAMETERS;

  int N = No_of_MODEL_PARAMETERS + No_of_IC + No_of_ERROR_PARAMETERS;
  
  F->Solution = (gsl_vector **)calloc(Realizations, sizeof(gsl_vector *));
  for(i = 0; i<Realizations; i++)
    F->Solution[i] = gsl_vector_alloc(N);

  F->Solution_Fitness = gsl_vector_alloc(Realizations);

  F->Solution_Order_Index = (int *)calloc(Realizations, sizeof(int) );

  F->No_of_SOLUTIONS = Realizations;

  F->TOTAL_No_of_Fitting_Parameters = N; 
}

void Parameter_Fitting_Initialization( Parameter_Fitting * F, int Realizations,
				       Observed_Data * Data, Parameter_Table * Table)
{
  int i;

  F->Data                  = Data;
  F->Space                 = Table->S;
  F->Table                 = Table;
  Table->Fitting_Data      = (void *)F;

  // Parameter Fitting and Parameter Table Structures point to each other!!! 
}

void Parameter_Fitting_Free(Parameter_Fitting * F)
{
  int i; 

  int Realizations = F->No_of_SOLUTIONS;
  
  for(i = 0; i<Realizations; i++)
    gsl_vector_free(F->Solution[i]);
  free(F->Solution);

  gsl_vector_free(F->Solution_Fitness);

  free(F->Solution_Order_Index);  
}

void Parametric_Configurations_into_Fitting_Structure_from_File (Parameter_Fitting * F,
								 char * File_Name )
{
  int i, j;
  int N;
  
  double ** Data = (double **)calloc( F->No_of_SOLUTIONS, sizeof(double *) );
  for( i=0; i < F->No_of_SOLUTIONS; i++ )
    Data[i] = (double *)calloc( F->TOTAL_No_of_Fitting_Parameters + 1, sizeof(double) );
  
  Reading_Model_Parameters_from_File(File_Name, Data, &N, F->TOTAL_No_of_Fitting_Parameters); 

  F->No_of_SOLUTIONS = N;
  
  for( i=0; i<N; i++ ) { 
    for( j=0; j < F->TOTAL_No_of_Fitting_Parameters; j++ )
      gsl_vector_set(F->Solution[i], j, Data[i][j]);  

    assert( j == F->TOTAL_No_of_Fitting_Parameters );
    
    gsl_vector_set(F->Solution_Fitness, i, Data[i][j]);
  }

  for( i=0; i < F->No_of_SOLUTIONS; i++ ) free(Data[i]);
  free(Data);
}

void Parametric_Configurations_from_Fitting_Structure_into_File (Parameter_Fitting * F,
								 char * File_Name,
								 int ORDERED_TRUE )
{
  int i, j;
  int key; 
  int N;
  double value; 
  gsl_vector * xx;     
  gsl_permutation * p;

  Parameter_Table * Table = F->Table;
  
  double ** Data = (double **)calloc( F->No_of_SOLUTIONS, sizeof(double *) );
  for( i=0; i < F->No_of_SOLUTIONS; i++ )
    Data[i] = (double *)calloc( F->TOTAL_No_of_Fitting_Parameters + 1, sizeof(double) );

  N = F->No_of_SOLUTIONS;

  if( ORDERED_TRUE == 1 ) {
    xx = gsl_vector_alloc(N);
    p  = gsl_permutation_alloc(N);

    for(i=0; i<N; i++) {
      value = gsl_vector_get(F->Solution_Fitness, i);
      gsl_vector_set(xx, i, value);
    }  
    gsl_sort_vector_index (p, xx);  
  }
    
  for( i=0; i<N; i++ ) {
    
    if ( ORDERED_TRUE == 1 ) key = gsl_permutation_get(p, i);
    else                     key = i;
    
    for( j=0; j < F->TOTAL_No_of_Fitting_Parameters; j++ ) {   
      Data[i][j] = gsl_vector_get(F->Solution[key], j);  
    }
    
    assert( j == F->TOTAL_No_of_Fitting_Parameters );
      
    Data[i][j] = gsl_vector_get(F->Solution_Fitness, key);
  }

  if( ORDERED_TRUE == 1 ) {
    gsl_vector_free(xx);
    gsl_permutation_free(p);
  }
   
  /* Creating Title Row */
  char ** Title_Parameters = (char **)calloc(F->TOTAL_No_of_Fitting_Parameters+1,
					     sizeof(char *) );
  for(i=0; i < F->TOTAL_No_of_Fitting_Parameters+1; i++)
    Title_Parameters[i] = (char *)calloc(100, sizeof(char));

  Creating_Title_Row (Table, Title_Parameters);
  
  Writing_Model_Parameters_into_File(File_Name, Title_Parameters, Data, N,
  				     F->TOTAL_No_of_Fitting_Parameters);

  for(i=0; i < j; i++) free(Title_Parameters[i]);
  free(Title_Parameters);

  for( i=0; i < F->No_of_SOLUTIONS; i++ ) free(Data[i]);
  free(Data);
}

void Accuracy_Calculation_from_Optimal_Parameter_Configuration(Parameter_Fitting * F,
							       int Input_Parameter_0,
							       double Parameter_True_Value_0,
							       double *Parameter_Estimated_Value_0,
							       double *Parameter_Accuracy_0,
							       int Input_Parameter_1,
							       double Parameter_True_Value_1,
							       double *Parameter_Estimated_Value_1,
							       double * Parameter_Accuracy_1)
{
  int i, j;
  int key; 
  int N;
  double value; 
  gsl_vector * xx;     
  gsl_permutation * p;
  
  Parameter_Table * Table = F->Table;
  
  double ** Data = (double **)calloc( F->No_of_SOLUTIONS, sizeof(double *) );
  for( i=0; i < F->No_of_SOLUTIONS; i++ )
    Data[i] = (double *)calloc( F->TOTAL_No_of_Fitting_Parameters + 1, sizeof(double) );

  N = F->No_of_SOLUTIONS;
  
  xx = gsl_vector_alloc(N);
  p  = gsl_permutation_alloc(N);
  
  for(i=0; i<N; i++) {
    value = gsl_vector_get(F->Solution_Fitness, i);
    gsl_vector_set(xx, i, value);
  }  
  gsl_sort_vector_index (p, xx);  

    
  for( i=0; i<N; i++ ) {
    
    key = gsl_permutation_get(p, i);
    
    for( j=0; j < F->TOTAL_No_of_Fitting_Parameters; j++ ) {   
      Data[i][j] = gsl_vector_get(F->Solution[key], j);  
    }
    
    assert( j == F->TOTAL_No_of_Fitting_Parameters );
      
    Data[i][j] = gsl_vector_get(F->Solution_Fitness, key);
  }

  gsl_vector_free(xx);
  gsl_permutation_free(p);

  /* Calculationg Accuracies */
  gsl_vector * x  = gsl_vector_alloc(F->Table->TOTAL_No_of_MODEL_PARAMETERS);
  for(i=0; i<F->Table->TOTAL_No_of_MODEL_PARAMETERS; i++)
    gsl_vector_set(x, i, Data[0][i]);
    
  Vector_Entries_into_Parameter_Table ( x, F->Table,
					Table->Index,
					F->Table->TOTAL_No_of_MODEL_PARAMETERS );

  * Parameter_Estimated_Value_0 = AssignStructValue_to_VectorEntry(Input_Parameter_0,
									F->Table);

  * Parameter_Estimated_Value_1 = AssignStructValue_to_VectorEntry(Input_Parameter_1,
									F->Table);

  * Parameter_Accuracy_0 = fabs( * Parameter_Estimated_Value_0 - Parameter_True_Value_0 ); 
  * Parameter_Accuracy_1 = fabs( * Parameter_Estimated_Value_1 - Parameter_True_Value_1 ); 

  /* Relative Error (with respect to true values) */
  * Parameter_Accuracy_0 /= Parameter_True_Value_0; 
  * Parameter_Accuracy_1 /= Parameter_True_Value_1; 
  
  printf("%s: Parameter True Value = %g\tParameter Estimated Value = %lf\tAccuracy = %lf\n",
	   Table->Symbol_Parameters[Input_Parameter_0],
	   Parameter_True_Value_0, * Parameter_Estimated_Value_0, * Parameter_Accuracy_0); 
  printf("%s: Parameter True Value = %g\tParameter Estimated Value = %lf\tAccuracy = %lf\n",
	   Table->Symbol_Parameters[Input_Parameter_1],
	   Parameter_True_Value_1, * Parameter_Estimated_Value_1, * Parameter_Accuracy_1); 

  gsl_vector_free( x ); 
  for( i=0; i < F->No_of_SOLUTIONS; i++ ) free(Data[i]);
  free(Data);
}
