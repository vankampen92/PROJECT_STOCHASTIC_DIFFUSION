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

#define EPSILON 1.0E-06

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
    #if defined STATIONARY_POINT_REPRESENTATION
    assert(P->No_of_CELLS == 1); /* Stationary Points work only for single patch systems */
    /* If Fixed Points have been calculated, comment out the following line:  */
    // Fixed_Points_All(Table, y_INI, y_INI, y_INI, EPSILON);  /* Calculating Lower/Single point... */
		//
    // Otherwise, if Fixed_Points_all is called in this function, you don't need the following two lines: */
    for( i = 0; i < Table->MODEL_STATE_VARIABLES; i++ )                    
         y_INI[i] = Table->Vector_Model_Variables_Time_0[i];
    #endif          
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

void Initial_Conditions_Global_Populations (Parameter_Table * Table, double * Vector_Model_Variables )
{
  /* Global populations across the patch system... checking for consistency */
  /* This function is calling structure at this level... 
     main() (main.c) 
          --->  M_O_D_E_L___S_T_O( &Table ); (MODEL_STO.c) 
                --->  S_T_O_C_H_A_S_T_I_C___T_I_M_E___D_Y_N_A_M_I_C_S (...); (Stochastic_Time_Dynamics.c)
                      ---> Initial_Conditions_Stochastic_Dynamics(...) (Initial_Conditions_Stochastic_Dynamics.c)
                      ---> Initial_Conditions_Global_Populations(...); (Initial_Conditions_Stochastic_Dynamics.c) 

     When called at this level, it yields global information about the strain and plasmid populations
     of the initial configuration of the system. 

     However, this call can be done at any point, since 'Vector_Model_Variables[]' specifies the 
     current configuration for the whole patch system at any given time (in particular, at Time 0, 
     but then it is updated for any further time).                    
  */
  int i, j, k, l, n;
  double N; 
  double * np = (double *)calloc(Table->No_of_RESOURCES, sizeof(double)); 

  assert(Table->TYPE_of_MODEL == 21 || Table->TYPE_of_MODEL == 20);

  /* Per Species Type (strain and profile in TYPE of MODEL 21) */
  for(k=0; k<Table->No_of_RESOURCES; k++) {
    Table->Focal_Resource = k;
    np[k] = Total_Population_Resource_Species (Vector_Model_Variables, Table); /* across cells */
    printf("Global Population (ID=%d) = %g\n", k, np[k]);
  }
  printf("\n");

  if(Table->TYPE_of_MODEL == 21) {
    /* Per Strain (suming over profiles) */
    j = 0;
    for(k=0; k<Table->No_of_STRAINS; k++) {
    
      N = 0.0; 
      for(i=0; i<Table->n[k]; i++) 
        N += np[j++];

      Table->Global_Strains_Population[k] = N;
      printf("Global Population Strain (over profiles) (Strain ID=%d) = %g\n", k, N); 
    }

    assert(j == Table->No_of_RESOURCES);

    /* Per Plasmid (suming over profiles and strain): 
        How many individual bacterial cells carry a given plasmid type? */
    for(k=0; k<Table->No_of_PLASMIDS; k++) {
    
      N = 0.0; l = 0; 
      for(j=0; j<Table->No_of_STRAINS; j++) {
        for(i=0; i<Table->n[j]; i++) {
          if( Table->Strain_Profiles[j][i][k] == 1 ) {
            N += np[l];
          }
          l++;
        }
      } 
      Table->Global_Plasmid_Population[k] = N;

      assert( l == Table->No_of_RESOURCES );  
    }
   
    printf("\n");
    /* Check: The same calculation, but explicitly using the patch system. Also, initialization:
              1. No of individual bacteria per Strain (regardless profile) for every cell in the patch systeem  
              2. No of individual bacteria carrying the same plasmid type for every cell in the patch system 
    */
    Community ** PATCH = Table->Patch_System;

    double * Total_Strain_Population = (double *)calloc(Table->No_of_STRAINS, sizeof(double));

    for(k=0; k<Table->No_of_CELLS; k++) {

      N = 0.0; l = 0; 
      for(j=0; j<Table->No_of_STRAINS; j++) {
        for(i=0; i<Table->n[j]; i++) {
          N += (double)PATCH[k]->Local_Strain_Population[l++]->n;
        }
        /* Local total number of individual cells per bacterial type (regardless profile) in the k-th Cell */
        PATCH[k]->Bacterial_Type_Population[j] = N;
        /* Global total number of individuals per bacterial type (across the system regardless profile) */
        Total_Strain_Population[j] += N; 
      }
    }

    for(j=0; j<Table->No_of_STRAINS; j++) {
      printf("Total Strain Population (across the system): Strain[%d] = %g\tStrain[%d] = %g\t", 
      j, Total_Strain_Population[j], 
      j, Table->Global_Strains_Population[j]);
    }
    free(Total_Strain_Population);

    double * Total_Plasmid_Population = (double *)calloc(Table->No_of_PLASMIDS, sizeof(double));
    for(k=0; k<Table->No_of_CELLS; k++) {
      
      for(n=0; n<Table->No_of_PLASMIDS; n++) { 

        l = 0; N = 0.0; 
        for(j=0; j<Table->No_of_STRAINS; j++) {
          for(i=0; i<Table->n[j]; i++) {
            if( Table->Strain_Profiles[j][i][n] == 1 ) {
              N += (double)PATCH[k]->Local_Strain_Population[l++]->n;
            }
          }
        }
        /* Local number of plasmid types (across all different bacterial types) */
        PATCH[k]->Plasmid_Type_Population[n] = N; 
        /* Global number of plasmid types (across the system for all different bacterial types) */    
        Total_Plasmid_Population[n] += N; 
      }
    }

    for(j=0; j<Table->No_of_PLASMIDS; j++) {
      printf("Total Plasmid Population (across the system): Plasmid[%d] = %g\tPlasmid[%d] = %g\t", 
      j, Total_Plasmid_Population[j], 
      j, Table->Global_Plasmid_Population[j]);
    }
    free(Total_Plasmid_Population);    
  }

  free(np);
  printf(" Initial Configuration has been checked. Checked passed.\n");
}

double Bacterial_Type_Population_per_Cell (Parameter_Table * Table, int i_Strain, int j_Patch)
{
  /* Local total number of individual cells per bacterial type or strain 'i_Strain 
     (regardless profile) in PATCH 'j_Patch' 
  */
  double N; 
  int i, l; 

  Community ** PATCH = Table->Patch_System;

  assert(i_Strain >= 0 && i_Strain < Table->No_of_STRAINS);
  
  l = 0;
  if(i_Strain > 0) 
    for(i=0; i<Table->n[i_Strain-1]; i++) 
      l += Table->n[i];
  
  N = 0.0;
  for(i=l; i<l+Table->n[i_Strain]; i++) {
          N += (double)PATCH[j_Patch]->Local_Strain_Population[i]->n;
  }

  return (N);
}

double Plasmid_Type_Population_per_Cell (Parameter_Table * Table, int i_Plasmid, int j_Patch)
{ 
  /* How prevalent is a plasmid in a given patch 'j_Patch'?
     Local number of plasmids (across all different bacterial types). 
     
     In particular, this functions responds to the questions: 
     How many bacterial individual cells carry the same plasmid 'i_Plasmid' 
     in the 'j_Patch' local patch, cell or local community?
  */

  double N; 
  int i, j, l; 

  Community ** PATCH = Table->Patch_System;

  l = 0; N = 0.0; 
  for(j=0; j<Table->No_of_STRAINS; j++) {
    for(i=0; i<Table->n[j]; i++) {
      if( Table->Strain_Profiles[j][i][i_Plasmid] == 1 ) {
        N += (double)PATCH[j_Patch]->Local_Strain_Population[l++]->n;
      }
    }
  }
  
  return (N);
}   