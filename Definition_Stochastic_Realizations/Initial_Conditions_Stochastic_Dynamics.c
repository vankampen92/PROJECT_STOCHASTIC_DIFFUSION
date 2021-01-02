#include <MODEL.h>

/*
   This function poses the PATCH system at the initial configuration and calculates
   the initial rate at which the total system may change from this initial configuration
   to any possible alternative one.
*/
#ifndef VERBOSE
#define STAT_BOOL_VERBOSE 0
#else
#define STAT_BOOL_VERBOSE 1
#endif

extern gsl_rng * r;

void Initial_Conditions_Stochastic_Dynamics( Parameter_Table * Table, double * y_INI )
{
  int i;

  Parameter_Model * P = Table->P;
  Community ** PATCH  = Table->Patch_System;

  if(Table->TYPE_of_INITIAL_CONDITION == 0) {
    Print_Press_Key (STAT_BOOL_VERBOSE, 0,
	  "Initial Conditions are defined from command line or by default values\n");
    Initial_Condition_from_Parameter_Table(Table, y_INI);
  }
  else if (Table->TYPE_of_INITIAL_CONDITION == 1) {
    Print_Press_Key (STAT_BOOL_VERBOSE, 0,
	  "Random Initial Conditions are defined within some boundary values\n");
    Random_Initial_Condition(Table, y_INI);
    Get_Initial_y_INI_Random_Vector_into_Integer_Values(Table, y_INI);
  }
  else if (Table->TYPE_of_INITIAL_CONDITION == 2) {
    Print_Press_Key (STAT_BOOL_VERBOSE, 0,
	  "Initial Conditions are defined as the fixed points of the system\n");
    // fixed_Points(P, y_INI, EPSILON);  /* Calculating Lower/Single point... */
				         /* (see fixed_Points.c)              */

    /* for( i = 0; i < Table->MODEL_STATE_VARIABLES; i++ )                    */
    /*    Table->Vector_Model_Variables_Stationarity[i] = y_INI[i] ;          */
  }
  else {
    printf(" Attention: Initial Condition Value is undefined:\n");
    printf(" Allows values are 0, 1, and 2, but TYPE_of_INITIAL_CONDITION = %d\n",
	   Table->TYPE_of_INITIAL_CONDITION );
    exit(0);
  }

  Patch_System_Initialization (PATCH, Table, y_INI);

  for (i=0; i<Table->MODEL_STATE_VARIABLES; i++ ) {
    Table->Vector_Model_Int_Variables[i]     = (int)y_INI[i];
    Table->Vector_Model_Variables_Time_0[i]  = (int)y_INI[i];
  }
}

void Get_Initial_y_INI_Random_Vector_into_Integer_Values(Parameter_Table * Table,
							 double * y_INI)
{
  int i,j,k,n;
  double value;

  int N = Table->INITIAL_TOTAL_POPULATION;
  int K = Table->MODEL_STATE_VARIABLES;
  int M = Table->No_of_CELLS;
  int S = Table->No_of_RESOURCES;

  printf("Initial No of Individuals per Species (input argument): %d\n", N);

  if (N < M) {
    double ** x = (double **)calloc(S, sizeof(double *));
    for(i=0; i<S; i++)
      x[i] = (double *)calloc(M, sizeof(double));
    gsl_vector * xx      = gsl_vector_alloc(M);
    gsl_permutation * p = gsl_permutation_alloc(M);

    for(i=0; i<S; i++)
      for(j=0; j<M; j++)
	x[i][j] = y_INI[i+j*S];

    for(i=0; i<S; i++) {

      for(j=0; j<M; j++)
	gsl_vector_set(xx, j, x[i][j]);

      gsl_sort_vector_index (p, xx);

      for(k=0; k < N; k++) {
	n = gsl_permutation_get(p, k);
	x[i][n] = 1.0;
      }
      for(k=N; k < M; k++){
	n = gsl_permutation_get(p, k);
	x[i][n] = 0.0;
      }
    }

    value = 0.0;
    for(i=0; i<S; i++)
      for(j=0; j<M; j++){
	y_INI[i+j*S] = x[i][j];
	value       += y_INI[i+j*S];
      }

    for(i=0; i<S; i++) free(x[i]);
    free(x);
    gsl_vector_free(xx);
    gsl_permutation_free(p);
  }
  else{
    printf("This random configuration does not make much sense!!!\n");
    printf("For a correct random initialization of the system, it is required\n");
    printf("that the number of cells or local populations is bigger than the\n");
    printf("the total initial population. Recall that the parameter INITIAL\n");
    printf("TOTAL POPULATION is a TOTAL INITIAL SIZE of the community\n");
    printf("(summing over the species in the system)\n");
    printf("FYI, INITIAL_TOTAL_POPULATION = %g\t No of CELLS = %d\n",
	   Table->INITIAL_TOTAL_POPULATION, Table->No_of_CELLS);
    printf("The program will exitn");
    exit(0);
  }

  printf("Comparing Total Initial Populations...\n");
  printf("Total Community Size (across species):  value = %g\n", value);
  printf("Total Population Size (per species, as input argument): %d\n", N);
  getchar();
  Table->INITIAL_TOTAL_POPULATION = value/(double)S;
}
