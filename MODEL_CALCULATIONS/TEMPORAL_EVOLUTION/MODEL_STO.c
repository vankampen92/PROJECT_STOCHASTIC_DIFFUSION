#include "../../Include/MODEL.h"

extern gsl_rng * r;

extern int TYPE_of_TIME_DEPENDENCE;

#if defined BINARY_TREE_OPTIMIZATION
    #define BINARY_TREE   
#endif
#if defined BINARY_TREE_SUPER_OPTIMIZATION
    #define BINARY_TREE   
#endif
#if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
    #define BINARY_TREE
#endif

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
  if (Table->No_of_CELLS == 1)
    if(Table->TYPE_of_MODEL == 12 || Table->TYPE_of_MODEL == 13 || 
       Table->TYPE_of_MODEL == 14 || Table->TYPE_of_MODEL == 16) {
      Initial_Condition_One_Single_Cell_into_Parameter_Table (Table,
                                  Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0,
                                  Table->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0);
    }
    else {
      Initial_Condition_One_Single_Cell_into_Parameter_Table (Table,
							                                                Table->INITIAL_TOTAL_POPULATION,
							                                                Table->INITIAL_TOTAL_POPULATION); 
                                                              /* -xN  (INITIAL_TOTAL_POPULATION) */
                                                              /* INITIAL_TOTAL_POPULATION::Initial_Conditions */
                                                              /* No_of_INDIVIDUALS::Parameter_Model           */                                  
    }
  else if(Table->No_of_CELLS <= 4)
    Initial_Condition_All_Patches_the_Same_into_Parameter_Table (Table,
								                                                 Table->INITIAL_TOTAL_POPULATION);
    /* All types and cells take the same initial total population */                                                          
  else 
    Initial_Condition_Centered_into_Parameter_Table (Table, 
                                                     Table->INITIAL_TOTAL_POPULATION); 
                                                     /* -xN [] */
  
  /* Initial condition: a stationary point (only when Table->No_of_CELLS is 1) */
  if (Table->TYPE_of_INITIAL_CONDITION == 2) {
    Print_Press_Key (1, 0, "Initial Conditions are defined as the fixed points of the 2D system\n");
    // Fixed Points should have been already calculated in a previous call 
    // to the function 'Fixed_Points_All();
    if (Table->No_of_CELLS == 1 && Table->TYPE_of_MODEL == 22) {
      for (i=0; i<MODEL_STATE_VARIABLES; i++) {
        Table->Vector_Model_Variables_Time_0[i] = Table->Vector_Model_Variables_Stationarity[i];
        
        // Table->Vector_Model_Variables_Time_0[i] = Table->Vector_Model_Variables_MultiStability[0][i];
        // Table->Vector_Model_Variables_Time_0[i] = Table->Vector_Model_Variables_MultiStability[0][i];
        // Table->Vector_Model_Variables_Time_0[i] = Table->Vector_Model_Variables_MultiStability[0][i];    
      }
    }
    else {
      printf("The system is not a single patch model (Table->No_of_CELLS = %d)\n", Table->No_of_CELLS);
      printf("The system is not a model of type 22 (Table->TYPE_of_MODEL = %d)\n", Table->TYPE_of_MODEL);
      printf("The program will safely exit\n");   
      exit(0);
    }
  }
  /* END ----------------------------------------------------------------------------
   */
  
  /* BEGIN : -------------------------------------------------------------------------
   * Stochastic Community Set Up
   */
  Community ** PATCH = (Community **)calloc( P->No_of_CELLS, sizeof(Community *) );
  Community_Allocation( PATCH, P ); 
  Community_Initialization (PATCH, P);
#if defined DIFFUSION_ECO_PLASMIDS
  Community_Plasmids_Initialization (PATCH, P);
  Community_Strains_Initialization (PATCH, P);  
#endif

  /* The Parameter Model structure also keeps the three memmory addresses pointing to 
   * the Patch System, the Time Control structure, and the CPG structure to plot   
   */
  Table->Patch_System = PATCH;

  #if defined BINARY_TREE_OPTIMIZATION
    /* See Community.c !!! */
    Community_Binary_Tree_Allocation (Table, Table->No_of_CELLS);   
  #endif
  #if defined BINARY_TREE_SUPER_OPTIMIZATION
    /* See Community.c !!! */
    Community_Binary_Tree_Allocation (Table, Table->TOTAL_GRAND_No_of_EVENTS);   
  #endif
  #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
    /* See Community.c !!! */
    Community_Priority_Queue_Tree_Allocation (Table, Table->TOTAL_GRAND_No_of_EVENTS); 
    // Print_Press_Key(1,1,"Printing out Tree before entering the generation of stochastic replicates");
    // printtree(Table->Treeroot);
  #endif
  #if defined BINARY_TREE
    P->Leaves            = Table->Leaves;
    P->No_of_LEAVES      = Table->No_of_LEAVES;                                 
    P->No_of_TREE_LEVELS = Table->No_of_TREE_LEVELS;
    /* Important assert: */
    assert(P->No_of_LEAVES == power_int(2, P->No_of_TREE_LEVELS));
   #endif 
  /* END ----------------------------------------------------------------------------
   */
  
#if defined CPGPLOT_REPRESENTATION  /* Initial Plotting Time evolution: just frames!!! */
  int SAME_PLOT = 0;
  // C_P_G___S_U_B___P_L_O_T_T_I_N_G___n___P_L_O_T_S( CPG->DEVICE_NUMBER,
  ////                                                SAME_PLOT, 0, Table );
#endif
  
  /* BEGIN: Main loop: a number of REALIZATIONS (stochastic temporal evolutions) is computed */
  printf("Entering Generation of Stochastic Realizations...\n");   
  Print_Press_Key(1,1,"Entering Generation of Stochastic Realizations...\n");
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
    int FROZEN_SYSTEM = S_T_O_C_H_A_S_T_I_C___T_I_M_E___D_Y_N_A_M_I_C_S ( n,
									                                                        Table, 
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
  /*         One file per each output variable            */  
  printf( " \nAs many as %d stochastic realizations have been successfully completed\n",
          Time->Realizations);
  printf( " Averages and Variances over the ensamble of realizations\n");
  printf( " will be calculated now...\n");
  Print_Press_Key(1,0,".");
  int DATA_POINTS = Time_Control_AVE_VAR_SAVE_VARIABLES( Table );
  printf(" Temporal series of %d (out of %d) data points\n", DATA_POINTS, I_Time);

#if defined CPGPLOT_REPRESENTATION
  // SAME_PLOT = 1; 
  // C_P_G___S_U_B___P_L_O_T_T_I_N_G___E_R_R_O_R___B_A_R ( Table, SAME_PLOT, 
  //							                                         DATA_POINTS, Table->T->time_DEF, 
  //							                                         Table->T->AVE, Table->T->VAR ); 
#endif
  /*   END : Averaging stochastic realizations            */  

  free( Table->Vector_Model_Variables );
  free( Table->Vector_Model_Variables_Time_0 ); 
  free( Table->Vector_Model_Int_Variables );
  free( Table->Vector_Model_Int_Variables_Time_0 );

#if defined STATIONARY_POINT_REPRESENTATION 
  // Fixed Points Calculations    
  /* De-allocating variables allocated in MODEL.c to calculate Fixed Points 
     that M_O_D_E_L___S_T_O( Parameter_Table * Table ) may require */
  free( Table->Vector_Model_Variables_MultiStability[0] );
  free( Table->Vector_Model_Variables_MultiStability[1] );
  free( Table->Vector_Model_Variables_MultiStability[2] );
  free( Table->Vector_Model_Variables_Stationarity );
#endif

#if defined DIFFUSION_ECO_PLASMIDS  
  Community_Strains_Free(PATCH, P);
#endif

  Community_Free(PATCH, P);
  free ( P ); 

  #if defined BINARY_TREE_OPTIMIZATION
    Binary_Tree_Free ( Table->Treeroot, Table->Leaves, Table->Parent, 
                       Table->No_of_CELLS ); 
  #endif
  #if defined BINARY_TREE_SUPER_OPTIMIZATION
    Binary_Tree_Free ( Table->Treeroot, Table->Leaves, Table->Parent, 
                       Table->TOTAL_GRAND_No_of_EVENTS ); 
  #endif
  #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
    Binary_Tree_Free ( Table->Treeroot, Table->Leaves, Table->Parent, 
                       Table->No_of_LEAVES );
    free(Table->Tree_Node_Index); /* Priority array of indexed pointers to all tree nodes */
  #endif 

  return(0);
}
