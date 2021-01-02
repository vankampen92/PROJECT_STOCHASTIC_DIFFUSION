#include <MODEL.h>

#if defined PARAMETRIC_CONFIGURATIONS
void Parameter_Space_Alloc( Parameter_Space * S, int No_of_PARAMETERS, int Realizations, int * d )
#else
void Parameter_Space_Alloc( Parameter_Space * S, int No_of_PARAMETERS, int * d )
#endif
{
  int i;

  S->P_MAX    = gsl_vector_alloc(No_of_PARAMETERS);
  S->P_min    = gsl_vector_alloc(No_of_PARAMETERS);
  S->Accuracy = gsl_vector_alloc(No_of_PARAMETERS);

  S->Parameter_Index =(int *)calloc(No_of_PARAMETERS, sizeof(int) );

  S->Parameter_MAX = (Parameter_Model *)calloc(1, sizeof(Parameter_Model) ) ;
  S->Parameter_min = (Parameter_Model *)calloc(1, sizeof(Parameter_Model) ) ;
  S->Parameter_Acc = (Parameter_Model *)calloc(1, sizeof(Parameter_Model) ) ;

  S->N = (int *)calloc(No_of_PARAMETERS, sizeof(int) );
  S->N_Par_Value = (double **)calloc(No_of_PARAMETERS, sizeof(double *) );
  for(i=0; i<No_of_PARAMETERS; i++)
    S->N_Par_Value[i] = (double *)calloc( d[i], sizeof(double) );

#if defined PARAMETRIC_CONFIGURATIONS
  int i;
  S->Parametric_Configuration = (double **)calloc( Realizations,
                                                   sizeof(double *) );
  for (i=0; i < Realizations; i++) {
    S->Parametric_Configuration[i] = (double *)calloc( No_of_PARAMETERS, sizeof(double) );
  }
#endif
}

void Parameter_Space_Initialization(Parameter_Space * S, int No_of_PARAMETERS,
				    double TOLERANCE, int MAX_No_of_ITERATIONS,
				    int * d, int * Index, double ** Range, double * Acc)
{
  int i, j, key;
  double lo_P, hi_P, ac_P;

  S->TOLERANCE            = TOLERANCE;
  S->MAX_No_of_ITERATIONS = MAX_No_of_ITERATIONS;
  S->No_of_PARAMETERS     =  No_of_PARAMETERS;

  for(i=0; i < S->No_of_PARAMETERS; i++)
    {  S->Parameter_Index[i] = Index[i]; S->N[i] = d[i]; }

  for(i=0; i < S->No_of_PARAMETERS; i++) {
    Vector_Entry_into_Parameter_Model( Range[0][i], S->Parameter_Index[i], S->Parameter_min );
    Vector_Entry_into_Parameter_Model( Range[1][i], S->Parameter_Index[i], S->Parameter_MAX );
    Vector_Entry_into_Parameter_Model( Acc[i],      S->Parameter_Index[i], S->Parameter_Acc );
  }

  for ( i=0; i < S->No_of_PARAMETERS; i++ ) {
    key = S->Parameter_Index[i];
    lo_P = Parameter_Model_into_Vector_Entry( key, S->Parameter_min );
    hi_P = Parameter_Model_into_Vector_Entry( key, S->Parameter_MAX );
    ac_P = Parameter_Model_into_Vector_Entry( key, S->Parameter_Acc );

    gsl_vector_set(S->P_min, i, lo_P);
    gsl_vector_set(S->P_MAX, i, hi_P);
    gsl_vector_set(S->Accuracy, i, ac_P);

    for(j=0; j<S->N[i]; j++)
      S->N_Par_Value[i][j] = lo_P + (double)j * (hi_P - lo_P)/(double)(S->N[i]-1);
  }
}

void Parameter_Space_IC_Initialization(Parameter_Space * S, int No_of_PARAMETERS,
				       int * d, int * Index, double ** Range, double * Acc)
{
  int i, j, key;
  double lo_P, hi_P, ac_P;

  S->No_of_PARAMETERS     =  No_of_PARAMETERS;

  for(i=0; i < S->No_of_PARAMETERS; i++)
    {  S->Parameter_Index[i] = Index[i]; S->N[i] = d[i]; }

  for(i=0; i < S->No_of_PARAMETERS; i++) {
    Vector_Entry_into_Model_Variable_Initial_Condition_Model( Range[0][i],
							      S->Parameter_Index[i],
							      S->Parameter_min );
    Vector_Entry_into_Model_Variable_Initial_Condition_Model( Range[1][i],
							      S->Parameter_Index[i],
							      S->Parameter_MAX );
    Vector_Entry_into_Model_Variable_Initial_Condition_Model ( Acc[i],
							       S->Parameter_Index[i],
							       S->Parameter_Acc );
  }

  for ( i=0; i < S->No_of_PARAMETERS; i++ ) {
    key = S->Parameter_Index[i];
    lo_P = Model_Variable_Initial_Condition_into_Vector_Entry_Model( key, S->Parameter_min );
    hi_P = Model_Variable_Initial_Condition_into_Vector_Entry_Model( key, S->Parameter_MAX );
    ac_P = Model_Variable_Initial_Condition_into_Vector_Entry_Model( key, S->Parameter_Acc );

    gsl_vector_set(S->P_min, i, lo_P);
    gsl_vector_set(S->P_MAX, i, hi_P);
    gsl_vector_set(S->Accuracy, i, ac_P);

    for(j=0; j<S->N[i]; j++)
      S->N_Par_Value[i][j] = lo_P + (double)j * (hi_P - lo_P)/(double)(S->N[i]-1);
  }
}

void Parameter_Space_Error_Initialization(Parameter_Space * S, int No_of_PARAMETERS,
				       int * d, int * Index, double ** Range, double * Acc)
{
  int i, j, key;
  double lo_P, hi_P, ac_P;

  S->No_of_PARAMETERS     =  No_of_PARAMETERS;

  for(i=0; i < S->No_of_PARAMETERS; i++)
    {  S->Parameter_Index[i] = Index[i]; S->N[i] = d[i]; }

  for(i=0; i < S->No_of_PARAMETERS; i++) {
    Vector_Entry_into_Error_Model_Model( Range[0][i],
     							      S->Parameter_Index[i],
							      S->Parameter_min );
    Vector_Entry_into_Error_Model_Model( Range[1][i],
							      S->Parameter_Index[i],
							      S->Parameter_MAX );
    Vector_Entry_into_Error_Model_Model ( Acc[i],
							       S->Parameter_Index[i],
							       S->Parameter_Acc );
  }

  for ( i=0; i < S->No_of_PARAMETERS; i++ ) {
    key = S->Parameter_Index[i];
    lo_P = Error_Model_into_Vector_Entry_Model( key, S->Parameter_min );
    hi_P = Error_Model_into_Vector_Entry_Model( key, S->Parameter_MAX );
    ac_P = Error_Model_into_Vector_Entry_Model( key, S->Parameter_Acc );

    gsl_vector_set(S->P_min, i, lo_P);
    gsl_vector_set(S->P_MAX, i, hi_P);
    gsl_vector_set(S->Accuracy, i, ac_P);

    for(j=0; j<S->N[i]; j++)
      S->N_Par_Value[i][j] = lo_P + (double)j * (hi_P - lo_P)/(double)(S->N[i]-1);
  }
}

void Parameter_Index_Checking_Ordering (int * A, int * D, int N,
					int * No_P0, int * No_P1, int * No_P2, int * No_P3,
					int * No_P4, int * No_P5, int * No_P6 )
{
  /* Input parameters:
     . A[i]: Definition of the parameter space through a list of integers
     containing the keys corresponding to model parameters.
     In general, the array I should be a permuation of the first N integers.
     If this is not true, the program aborts.
     . N: Length of arrays A[] and D[]: Maximum number of model parameters.
     . D: Parameter Space Discretization: No of discrete points in each direction of
          newly redefined (sub)parameter space.

     Output parameters:
     . No_C: Discretization in the first (key = 0) parameter
     . No_E: Discretization in the second (key = 1) parameter
     . No_D: Discretization in the third (key = 2) parameter
     . No_P: Discretization in the fourth (key = 3) parameter
  */
  int i,j,n,m;
  double * Index = (double *)calloc( N, sizeof(double) );
  size_t    * p     = (size_t   *)calloc( N, sizeof(size_t)    );

  m = 0;
  for(i=0; i<N; i++){
    n = 0;
    for (j=0; j<N; j++) if(A[j] == i) n++;
    if (n == 1) m++;

    Index[i] = (double)A[i];
  }
  assert( m == N );

  /* gsl_sort_index (...):
       This function indirectly sorts the n elements of the array 'Index' with stride 'double'
       into ascending order, storing the resulting permutation in p. The array p must be
       allocated with a sufficient length to store the n elements of the permutation. The
       elements of p give the index of the array element which would have been stored in that
       position if the array had been sorted in place. The array Index is not changed.
  */
    gsl_sort_index (p, Index, sizeof(double), N);

    (* No_P0) = D[p[0]];
    (* No_P1) = D[p[1]];
    (* No_P2) = D[p[2]];
    (* No_P3) = D[p[3]];
    (* No_P4) = D[p[4]];
    (* No_P5) = D[p[5]];
    (* No_P6) = D[p[6]];

    free(Index); free(p);
}

void Parameter_Space_Free(Parameter_Space * S, int No_of_PARAMETERS)
{
  int i;

  gsl_vector_free(S->P_MAX);
  gsl_vector_free(S->P_min);
  gsl_vector_free(S->Accuracy);
  free (S->Parameter_Index);
  free (S->Parameter_MAX);
  free (S->Parameter_min);
  free (S->Parameter_Acc);

  for(i=0; i<No_of_PARAMETERS; i++) free( S->N_Par_Value[i] );
  free( S->N_Par_Value );

  free( S->N );

  //free( S );
}

void Parameter_Space_Write_Min_Max_Error_Values (Parameter_Space * S, Parameter_Table * Table )
{
  int i, key;
  double lo_P, hi_P;
  int No_of_PARAMETERS  = S->No_of_PARAMETERS;

  printf("Parameter Space (Initial Conditions):\nSymbol\t Minimum Value\t Maximum Value\n");
  for(i = 0; i<No_of_PARAMETERS; i++){
    key = S->Parameter_Index[i];
    /* lo_P = Model_Variable_Initial_Condition_into_Vector_Entry_Model( key, S->Parameter_min ); */
    /* hi_P = Model_Variable_Initial_Condition_into_Vector_Entry_Model( key, S->Parameter_MAX ); */
    lo_P = gsl_vector_get(S->P_min, i);
    hi_P = gsl_vector_get(S->P_MAX, i);
    if(key < Table->OUTPUT_VARIABLES_GENUINE)
    printf("%s\t\t%g\t%g\n", Table->Output_Variable_Symbol[key], lo_P, hi_P);
    else
    printf("%s\t\t%g\t%g\n", "Error Parameter", lo_P, hi_P);
  }
}

void Parameter_Space_Write_Min_Max_IC_Values (Parameter_Space * S, Parameter_Table * Table )
{
  int i, key;
  double lo_P, hi_P;
  int No_of_PARAMETERS  = S->No_of_PARAMETERS;

  printf("Parameter Space (Initial Conditions):\nSymbol\t Minimum Value\t Maximum Value\n");
  for(i = 0; i<No_of_PARAMETERS; i++){
    key = S->Parameter_Index[i];
    /* lo_P = Model_Variable_Initial_Condition_into_Vector_Entry_Model( key, S->Parameter_min ); */
    /* hi_P = Model_Variable_Initial_Condition_into_Vector_Entry_Model( key, S->Parameter_MAX ); */
    lo_P = gsl_vector_get(S->P_min, i);
    hi_P = gsl_vector_get(S->P_MAX, i);

    printf("%s\t\t%g\t%g\n", Table->Model_Variable_Symbol[key], lo_P, hi_P);
  }
}

void Parameter_Space_Write_Min_Max_Values (Parameter_Space * S, Parameter_Table * Table )
{
  int i, key;
  double lo_P, hi_P;
  int No_of_PARAMETERS  = S->No_of_PARAMETERS;

  printf("Parameter Space:\nSymbol\t Minimum Value\t Maximum Value\n");
  for(i = 0; i<No_of_PARAMETERS; i++){
    key = S->Parameter_Index[i];
    /* lo_P = Parameter_Model_into_Vector_Entry( key, S->Parameter_min ); */
    /* hi_P = Parameter_Model_into_Vector_Entry( key, S->Parameter_MAX ); */

    lo_P = gsl_vector_get(S->P_min, i);
    hi_P = gsl_vector_get(S->P_MAX, i);

    printf("%s\t\t%g\t%g\n", Table->Symbol_Parameters[key], lo_P, hi_P);
  }
}
