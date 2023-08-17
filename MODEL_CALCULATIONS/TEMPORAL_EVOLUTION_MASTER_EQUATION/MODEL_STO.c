#include "../../Include/MODEL.h"

extern gsl_rng * r;

extern int TYPE_of_TIME_DEPENDENCE;

int M_O_D_E_L___S_T_O( Parameter_Table * Table )
{
  int i,j,k, n;
  int I_Time, no_Patch;
  int Bad_Times;
  double t; 
  Time_Control * Time;
  Time_Dependence_Control * TDC; 

  Time = Table->T;
  
  TDC  = Table->TDC; 
  
  Parameter_Model * P = (Parameter_Model *)malloc( 1 * sizeof(Parameter_Model) );
  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N (Table, P);
  Table->P  = P;
  P->CPG    = Table->CPG_STO; 
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
    if(Table->TYPE_of_MODEL == 12 || Table->TYPE_of_MODEL == 13 || Table->TYPE_of_MODEL == 14 || Table->TYPE_of_MODEL == 16) {
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
  /* The Parameter Table (and Parameter Model) structures also keep the three memmory addresses 
   * pointing to the Patch System (Parameter_Model), the Time Control structure, the CPG plotting 
   * structure, and the two structures controling the binary tree of 'treenode' type: 'Treeroot' and 'Leaves' 
   * (also stored in Parameter Model).    
   */
  Table->Patch_System = PATCH;
  
  #if defined BINARY_TREE_OPTIMIZATION
    Community_Binary_Tree_Allocation (Table);   /* See Community.c !!!  */
    /* Allocation of the Binary Tree in Memmory (with zero values) !!! */
    P->Leaves   = Table->Leaves; 
  #endif
  /* END ----------------------------------------------------------------------------
  */
  
#if defined CPGPLOT_REPRESENTATION  /* Initial Plotting Time evolution: just frames!!! */
  int SAME_PLOT = 0;
  // C_P_G___S_U_B___P_L_O_T_T_I_N_G___n___P_L_O_T_S( CPG->DEVICE_NUMBER,
  //                                                  SAME_PLOT, 0, Table );
#endif
  
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
       . n, lable of current realization 
       . P, a comprehensive model parameter table (see definition in MODEL.h)
       . Bad_Times is a measure of the performance of the sampling frequency. 
         If Bad_Times is high, interval times should be choosen smaller 
    */
    int FROZEN_SYSTEM = S_T_O_C_H_A_S_T_I_C___T_I_M_E___D_Y_N_A_M_I_C_S ( n, Table, 
                                                                          &Bad_Times );
    /* End of the i-th STOCHASTIC REALIZATIONS */
    printf("Realization: %d of a total of %d\n", n+1, Table->T->Realizations);
    printf("Time failed in %d occasions out of %d time steps\n", Bad_Times, I_Time);
    printf("If the number of failed times is too big, EPSILON might be too small!\n");
    printf("Try to choose a larger EPSILON [Current value: -tE %g]\n", Table->T->EPSILON);

    /* Only selection those stochastic realizations according to certain criterion */
    if(FROZEN_SYSTEM == 0) n++;  
  }
  /* END: End of STOCHASTIC REALIZATIONS */

  /* BEGIN : Averaging and saving stochastic realizations */
  printf( " \nAs many as %d stochastic realizations have been successfully completed\n",
          Time->Realizations);
  printf( " Averages and Variances over the ensamble of realizations\n");
  printf( " will be calculated now...\n");
  Print_Press_Key(1,0,".");
  int DATA_POINTS = Time_Control_AVE_VAR_SAVE_VARIABLES( Table );
  printf(" Temporal series of %d (out of %d) data points\n",
	 DATA_POINTS, I_Time);

#if defined CPGPLOT_REPRESENTATION
  SAME_PLOT = 1; 
  C_P_G___S_U_B___P_L_O_T_T_I_N_G___E_R_R_O_R___B_A_R ( Table, SAME_PLOT, 
							                                          DATA_POINTS, Table->T->time_DEF, 
							                                          Table->T->AVE, Table->T->VAR ); 
#endif
  /*   END : Averaging stochastic realizations -------------------------*/  
  
  free( Table->Vector_Model_Variables );
  free( Table->Vector_Model_Variables_Time_0 ); 
  free( Table->Vector_Model_Int_Variables );
  free( Table->Vector_Model_Int_Variables_Time_0 );
  
  Community_Free(PATCH, P);
  free ( P ); 

  #if defined BINARY_TREE_OPTIMIZATION
    Community_Binary_Tree_Free (Table);   /* See Community.c !!!  */
  #endif
    
  return(0);
}

