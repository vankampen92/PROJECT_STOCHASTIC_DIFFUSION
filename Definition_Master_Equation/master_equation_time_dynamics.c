#include <MODEL.h>
/*
   This driver produces the temporal evolution of a 
   master equation system. Sampling times are defined in 
   Time structure when this data structure is setup from 
   scratch at the main program.
   
   PGPLOTing is also possible at running time.

   This function calls the gsl-based ODE driver:

   int master_equation_driver( Parameter_Table * Table, int i, double * Time_Current );

   which is a generic common function which is always called
   for any implemented particular model
*/

int master_equation_time_dynamics( Parameter_Table * Table )
{
  /* This function performs a numerical integration of the master equation for a system 
     advancing from a time to the next time as stored in Time->Time_Vector[], and saves a 
     file corresponding to this numerical integration.

     This function makes two essential calls:

     1. Initial_Condition_Master_Equation (...), which sets up initial conditions. 

     2. master_equation_driver(...), which performs the actual generic numerical integration

     In addition, it depends on the following functions:
     1. Time_Dependence_Apply (...), which, in turn, calls Time_Dependence_Function(...),
     which performs the update of parameter table depending on the value of the input parameter
     TYPE_of_TIME_DEPENDENCE

     2. definition_OutPut_Variables(...)
  */
  int i;
  int State;
  // FILE *FP; char file[50];
  int j, k, kk;
  int TIMES;
  Time_Control * Time;
  double Time_Current, value;

  Time         = Table->T;
  TIMES        = Table->T->I_Time;

  // file[0]='\0';
  // fitxer(file, "ME_output_", 0, ".dat"); FP = fopen(file, "w");

  /* BEGIN : Initial Conditions * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#if defined VERBOSE
  printf(" Before Initial Condition Master Equation (...)\n");
#endif
  
  Initial_Condition_Master_Equation( Table, Table->MEq->Probability_Distribution );

#if defined VERBOSE
  printf(" After Initial Condition Master Equation (...)\n");
#endif
  
  for (k=0; k < Table->MEq->No_of_CONFIGURATIONAL_STATES; k++) {
    Table->MEq->Probability_Distribution_Time_0[k] = Table->MEq->Probability_Distribution[k];
  }
  
  Time_Current = Time->Time_Vector[0];
  if (Time->TYPE_of_TIME_DEPENDENCE > 0) Time_Dependence_Apply( Table, Time_Current );
  /*   END : Initial Conditions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#if defined VERBOSE
  printf(" Initiating Numerical Integration Master Equation \n");
#endif
  Marginal_Probability_Calculation ( Table );           /* At time zero              */
  Marginal_Probability_Averages_Calculation ( Table );  /* At time zero: <n> and <m> */
  Print_Probability_Distribution ( Table );
  Print_Marginal_Averages (Time_Current, Table);

#if defined CPGPLOT_REPRESENTATION
  int SAME = 0;
  j = 0;
  for ( i=0; i< Table->MEq->n_DIMENSION; i++ ) {  
    SAME = 0; 
    C_P_G___M_A_R_G_I_N_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Table, j, i, 
                                                        Time_Current,
							                                          SAME);
    SAME = 1; 
    C_P_G___T_H_E_O_R_E_T_I_C_A_L___M_A_R_G_I_N_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Table, j, i, 
                                                                                Time_Current,
	  						                                                                SAME );
    /* BEGIN: Saving Marginal for current time */
    Saving_Marginal_Distribution(Table, j, i, Time_Current);
    /*   END: -------------------------------- */ 
  }
  #if defined DIFFUSION_BD_2D    
    Saving_Marginal_Distribution_Triplets(Table, j, Time_Current);
  #endif 
  #if defined VERBOSE
    Print_Press_Key(0, 1,".");
  #endif
#endif

  for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
    kk = Table->OUTPUT_VARIABLE_INDEX[k];
    value = definition_OutPut_Variables( kk,
					                               Table->MEq->Vector_Model_Variables, 
                                         Time->Time_Vector[0],
					                               Table );
    Table->Matrix_Output_Variables[k][0] = value;
  }
  
  /* Main Loop: Advancing Time accoding to a Time (predefined) Vector */
  SAME = 0;
  for( j = 1; j < TIMES; j++ ) {
    /* This loop advances the system sequentially from
       intitial time 0 to 1st time , ...,  from time (j-1) to j, and so on.
       Note: When the system is frozen (FROZEN_SYSTEM = 1), then
       this loop does not advance the system any more
    */
    if (Table->T->TYPE_of_TIME_DEPENDENCE > 0) Time_Dependence_Apply( Table, Time_Current );
    
    /*-------------------------------------------------------------------*/
    /* B E G I N :
     *  CORE POINT HERE: Numerical Integration of the master equation up to the next time
     */
      State = master_equation_driver( Table, j, &Time_Current );
    /* Normalization: */
      Normalization_Master_Equation(Table);
    /* E N D : ------------------------------------------------------*/

    if (State != GSL_SUCCESS) break;

  #if defined CPGPLOT_REPRESENTATION
    /* This should be only activated in case we want to animate ODE time evolution by
       representing the solution as time progresses                                       
    */
    Marginal_Probability_Calculation ( Table );
    Marginal_Probability_Averages_Calculation ( Table );  /* At time zero */
    for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
      kk = Table->OUTPUT_VARIABLE_INDEX[k];
      value = definition_OutPut_Variables( kk,
					                                 Table->MEq->Vector_Model_Variables, 
                                           Time_Current, 
					                                 Table );
      Table->Matrix_Output_Variables[k][j] = value;
    }

    for ( i=0; i< Table->MEq->n_DIMENSION; i++ ) {
      SAME = 0; 
      C_P_G___M_A_R_G_I_N_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Table, j, i, 
                                                          Time_Current, 
                                                          SAME );
      SAME = 1; 
      C_P_G___T_H_E_O_R_E_T_I_C_A_L___M_A_R_G_I_N_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Table, j, i,  
			 				                                                                    Time_Current,
			 				                                                                    SAME );
      /* BEGIN: Saving Marginal for current time */
      Saving_Marginal_Distribution(Table, j, i, Time_Current);
      /*   END: -------------------------------- */
    }
    #if defined DIFFUSION_BD_2D    
      Saving_Marginal_Distribution_Triplets(Table, j, Time_Current);
    #endif   
  #endif 

    Print_Press_Key(0, 1,".");
  }/* ------> go further to the next time step */

#if defined VERBOSE
  printf(" Numerical Integration of the Master Equation ended!!!\n");
#endif

  // fclose(FP);

#if defined VERBOSE
  if (State != GSL_SUCCESS ) {
    printf(" Numerical integration failed at the %dth time\n", j);
  }
#endif
  
  return(State);
}
