#include <MODEL.h>

#define EPSILON 1.0E-06

extern gsl_rng * r;

int M_O_D_E_L( Parameter_Table * Table )
{
  double x;

  int i,k, n;
  int I_Time, no_Patch;
  int Bad_Times;
  Parameter_Model * P;
  Parameter_Fitting * F = (Parameter_Fitting *)Table->Fitting_Data;
  double ** Data        =  F->Data->N;

  I_Time    = Table->T->I_Time;

  P = (Parameter_Model *)malloc( 1 * sizeof(Parameter_Model) );
  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N (Table, P);
  Table->P  = P;
#if defined VERBOSE  
  printf(" Parameter_Model structure has been correctly allocated and initiated\n");
#endif
  
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  int MODEL_STATE_VARIABLES = K+1;
  Table->MODEL_STATE_VARIABLES = MODEL_STATE_VARIABLES;
  Table->Vector_Model_Variables_Time_0 = (double *)calloc( MODEL_STATE_VARIABLES, sizeof(double)); 
  Table->Vector_Model_Variables = (double *)calloc( MODEL_STATE_VARIABLES, sizeof(double) );
  
  /* BEGIN : -------------------------------------------------------------------------
   * Definition Initial Condition (initializing 'Table->Vector_Model_Variables_Time_0' vector):
   */
  if(Table->No_of_CELLS > 4)
    Initial_Condition_Centered_into_Parameter_Table (Table, Table->INITIAL_TOTAL_POPULATION);
  else if (Table->No_of_CELLS == 1)
    /* Initial Condition from Empirical Data File (stored now in a data structure, a member 
       of Parameter_Fitting structure, F 
    */
    /* The system has only one cell (no spatial structure) */
    Initial_Condition_One_Single_Cell_into_Parameter_Table (Table,
							    Data[0][0], 
							    Data[1][0]);
  else 
    Initial_Condition_All_Patches_the_Same_into_Parameter_Table (Table,
								 Table->INITIAL_TOTAL_POPULATION);
  /* END ----------------------------------------------------------------------------
   */
  /* BEGIN : -------------------------------------------------------------------------
   * Community Set Up
   */
  Community ** PATCH = (Community **)malloc( P->No_of_CELLS * sizeof(Community *) );
  Community_Allocation( PATCH, P ); 
  Community_Initialization (PATCH, P);
  /* The Parameter Model structure also keeps the three memmory addresses pointing to 
   *   the Patch System, the Time Control structure, and the CPG structure to plot   
   */
  Table->Patch_System = PATCH;
  /* END ----------------------------------------------------------------------------
   */			  							   
#if defined VERBOSE  
  printf(" Entering deterministic dynamics. Parameter time dependencies will be\n");
  printf(" de-activated if -t4 0 (TYPE_of_TIME_DEPENDENCE = 0).\n");
#endif
  D_E_T_E_R_M_I_N_I_S_T_I_C___T_I_M_E___D_Y_N_A_M_I_C_S( Table ) ;

#if defined CPGPLOT_REPRESENTATION
  //  Parameter Table dependent costumized plotting is defined
  // in ~/CPGPLOT/CPGPLOT_Parameter_Table/ directory
  /// int TIMES           = Table->T->I_Time;
  // int Input_Parameter = 0; /* The value of this model parameter appears in the title */
  // C_P_G___S_U_B___P_L_O_T_T_I_N_G ( Table, TIMES, Table->CPG->x_Time, Table->CPG->y_Time );
  // C_P_G___S_U_B___P_L_O_T_T_I_N_G___C_U_S_T_O_M_I_Z_E_D___T_I_T_L_E ( Table,
  //  								      TIMES,
  //  								      Table->CPG->x_Time,
  //  								      Table->CPG->y_Time,
  //  								      Input_Parameter );
#endif
  
  free( Table->Vector_Model_Variables_Time_0);
  free( Table->Vector_Model_Variables );

  Community_Free(PATCH, P);
  free ( P );
  
#if defined VERBOSE
  printf(" Deterministic dynamics successfully completed!!!\n");
#endif

  return(0);
}
