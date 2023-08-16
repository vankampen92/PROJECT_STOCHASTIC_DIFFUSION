#include "../../Include/MODEL.h"

extern gsl_rng * r;

extern int TYPE_of_TIME_DEPENDENCE;

int M_O_D_E_L___S_T_O( Parameter_Table * Table )
{
  /* Notice that the version of this function here is devoted to generate 
     stochastic replciates in a fully numerical way. It does not represent 
     anything on the fly. It only generates stochastic temporal replicates 
     of the dynamics specified by the MODEL environment variable and 
     for the specified output variables. 
  */
  int i,j,k, n;
  int I_Time, no_Patch;
  int Bad_Times;
  double t; 
  Time_Control * Time           = Table->T;
  Time_Dependence_Control * TDC = Table->TDC; 

  Parameter_Model * P = (Parameter_Model *)malloc( 1 * sizeof(Parameter_Model) );
  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N (Table, P);
  Table->P  = P;
  printf(" Parameter_Model structure has been correctly allocated and initiated\n");
  I_Time    = Table->T->I_Time;

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */  
  int MODEL_STATE_VARIABLES = Table->MODEL_STATE_VARIABLES;
  Table->Vector_Model_Variables = (double *)calloc( MODEL_STATE_VARIABLES, sizeof(double) );
  Table->Vector_Model_Variables_Time_0 = (double *)calloc( MODEL_STATE_VARIABLES, sizeof(double));
  Table->Vector_Model_Int_Variables = (int *)calloc( MODEL_STATE_VARIABLES, sizeof(int) );  
  Table->Vector_Model_Int_Variables_Time_0 = (int *)calloc( MODEL_STATE_VARIABLES, sizeof(int) );
  P->Vector_Model_Variables            = Table->Vector_Model_Variables;
  P->Vector_Model_Variables_Time_0     = Table->Vector_Model_Variables;
  P->Vector_Model_Int_Variables        = Table->Vector_Model_Int_Variables;    
  P->Vector_Model_Int_Variables_Time_0 = Table->Vector_Model_Int_Variables_Time_0;

  /* BEGIN : -------------------------------------------------------------------------
   * Definition Initial Condition (initializing 'Table->Vector_Model_Variables_Time_0' vector): 
   */
  if (Table->No_of_CELLS > 4)
    Initial_Condition_Centered_into_Parameter_Table (Table, Table->INITIAL_TOTAL_POPULATION);
  else if (Table->No_of_CELLS == 1)
    if(Table->TYPE_of_MODEL == 12 || Table->TYPE_of_MODEL == 13 || Table->TYPE_of_MODEL == 14 ) {
      Initial_Condition_One_Single_Cell_into_Parameter_Table (Table,
						   Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0,
						   Table->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0);
    }
    else {
      Initial_Condition_One_Single_Cell_into_Parameter_Table (Table,
							      Table->INITIAL_TOTAL_POPULATION,
							      Table->INITIAL_TOTAL_POPULATION);
    }
  else 
    Initial_Condition_All_Patches_the_Same_into_Parameter_Table (Table,
								 Table->INITIAL_TOTAL_POPULATION);
  
  for(i=0; i<Table->MODEL_STATE_VARIABLES; i++)
    Table->Vector_Model_Int_Variables_Time_0[i] = (int)Table->Vector_Model_Variables_Time_0[i];
  /* END ----------------------------------------------------------------------------
   */

  /* BEGIN : -------------------------------------------------------------------------
   * Stochastic Community Set Up
   */
  Community ** PATCH = (Community **)malloc( P->No_of_CELLS * sizeof(Community *) );
  Community_Allocation( PATCH, P ); 
  Community_Initialization (PATCH, P);
  /* The Parameter Model structure also keeps the three memmory addresses pointing to 
   * the Patch System, the Time Control structure, and the CPG structure to plot   
   */
  Table->Patch_System = PATCH;
  /* END ----------------------------------------------------------------------------
   */
    
  /* BEGIN: Main loop: a number of REALIZATIONS (stochastic temporal evolutions) is computed */
  printf("Entering Generation of Stochastic Realizations...\n");   Print_Press_Key(1,0,".");
  n = 0;
  while (n < Table->T->Realizations){
    // Notice that TDC has not been initialized when TYPE_of_TIME_DEPENDENCE = 0
    // This is the reason we need an 'extern int TYPE_of_TIME_DEPENDENCE' above!!!
    if (TYPE_of_TIME_DEPENDENCE == 1) {
      if(TDC->TYPE_2_PARAMETERS > 0) {
	for(i = 0; i < TDC->TYPE_2_PARAMETERS; i++) {
	  k = i + TDC->TYPE_0_PARAMETERS + TDC->TYPE_1_PARAMETERS;
	  for(j = 0; j<TDC->No_of_TIMES; j++) {
	    t = Time->Time_Vector[j];
	    TDC->Dependent_Parameter[k][j] =Time_Dependence_Resolve(Table,
								    TDC->Index_Dependent_Parameters[k],
								    TDC->Forcing_Pattern_Parameters[k], t);
	    
	  }
	}
      }
    }

    // GSL_Init_Random_Seed(r); /* According to Computer time */
    /* Input variables: 
       . i, lable of current realization 
       . P, a comprehensive model parameter table (see definition in MODEL.h)
       . Bad_Times is a measure of the performance of the sampling frequency. 
         If Bad_Times is high, interval times should be choosen smaller 
    */
    int FROZEN_SYSTEM = Stochastic_Time_Dynamics_Numerical ( n,
							     Table, &Bad_Times );
    
    /* End of the i-th STOCHASTIC REALIZATIONS */
    printf("Realization: %d of a total of %d\n", n+1, Table->T->Realizations);
    printf("Time failed in %d occasions out of %d time steps\n", Bad_Times, I_Time);
    printf("If the number of failed times is too big, EPSILON might be too small!\n");
    printf("Try to choose a larger EPSILON [Current value: -tE %g]\n", Table->T->EPSILON);

    /* Only selection those stochastic realizations according to certain criterion */
    if(FROZEN_SYSTEM == 0) n++;  

  }
  /* END: End of STOCHASTIC REALIZATIONS */

  printf( " \nAs many as %d stochastic realizations have been successfully completed\n",
          Time->Realizations);
  
  free( Table->Vector_Model_Variables );
  free( Table->Vector_Model_Variables_Time_0 ); 
  free( Table->Vector_Model_Int_Variables );
  free( Table->Vector_Model_Int_Variables_Time_0 );
  
  Community_Free(PATCH, P);
  free ( P ); 
  
  return(0);
}


