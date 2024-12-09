#include <MODEL.h>

#define EPSILON (1.0e-8)           /* Accuracy of the iterative process to find 
				      the stationary solution */

#ifndef VERBOSE
#define STAT_BOOL_VERBOSE 0
#else
#define STAT_BOOL_VERBOSE 1
#endif

extern gsl_rng * r; /* Global generator defined in main.c */

void Initial_Conditions_Numerical_Integration( Parameter_Table * Table, double * y_INI )
{
  int i; 
 
  if(Table->TYPE_of_INITIAL_CONDITION == 0) {
     
	  printf("Initial Conditions: alls species a the lowest abundance(see Initial_Conditions_Lower_Bound(...))\n");
    getchar(); 
    Initial_Condition_Lower_Bound( Table, y_INI );
  }
  else if (Table->TYPE_of_INITIAL_CONDITION == 1) {
  
	  printf("Initial Conditions are random between 0.0 and system size (or carring capacity K)\n");
    getchar(); 
    Random_Initial_Condition( Table, y_INI );
  }
  else if (Table->TYPE_of_INITIAL_CONDITION == 2) {
    
    printf("Initial Conditions are defined as the fixed points of the system\n");
    printf(" Not yet implemented...\n");
    exit(0); 
    // Fixed Points should have been already calculated in a previous call 
    // to the function 'Fixed_Points_All();
    // for( i = 0; i < Table->MODEL_STATE_VARIABLES; i++ ) 
    //   y_INI[i] = Table->Vector_Model_Variables_Stationarity[i];
  }
  else {
    printf(" Attention: Initial Condition Value is undefined:\n");
    printf(" Allows values are 0, and 1, but TYPE_of_INITIAL_CONDITION = %d\n",
	         Table->TYPE_of_INITIAL_CONDITION );
    exit(0); 
  }    
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
}

void Initial_Condition_Lower_Bound( Parameter_Table * Table, double * y_INI )							 
{
  int i; 

  for(i = 0; i<Table->No_of_RESOURCES; i++)
    y_INI[i] = 0.001 * (double)Table->K_R; 
}  

void Random_Initial_Condition( Parameter_Table * Table, double * y_INI )							 
{
  int i; 
  
  for(i = 0; i<Table->No_of_RESOURCES; i++)
    y_INI[i] = (double)Table->K_R * gsl_rng_uniform_pos(r); 
}   
  
