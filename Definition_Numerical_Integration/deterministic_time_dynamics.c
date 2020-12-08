#include <MODEL.h>
/*
   This driver produces the temporal evolution of a ODE
   system. Sampling times are defined in Time structure
   when this data structure is setup from scratch
   at the main program.
   PGPLOTing is also possible at running time.

   This function calls the gsl-based ODE driver:

   int numerical_Integration_Driver( Parameter_Table * Table, int i, double * Time_Current );

   which is a generic common function which is always called
   for any implemented particular model
*/
// #define VERBOSE

int D_E_T_E_R_M_I_N_I_S_T_I_C___T_I_M_E___D_Y_N_A_M_I_C_S( Parameter_Table * Table )
{
  /* This function performs a numerical integration of a the system avancing from a time 
     to the next as stored in Time->Time_Vector[], and save a file corresponding to this 
     numerical integration.

     This function makes two essential calls:

     1. Initial_Conditions_Numerical_Integration (...),
     which sets up initial conditions. In turn, this function, depending on the input parameter
     TYPE_of_INITIAL_CONDITION, sets up initial consitions in one way or another, by calling
     the functions:
     . Initial_Condition_from_Parameter_Table(P, y_INI);
     . Random_Initial_Condition(P, y_INI);
     . fixed_Points(P, y_INI, EPSILON);

     2. numerical_Integration_Driver(...), which performs the actual numerical integration

     In addition, it depends on the following functions:
     1. Time_Dependence_Apply (...), which, in turn, calls Time_Dependence_Function(...),
     which performs the update of parameter table depending on the value of the input parameter
     TYPE_of_TIME_DEPENDENCE

     2. definition_OutPut_Variables(...)
  */
  int i; int State;
  // FILE *FP; char file[50];
  int j, k, kk;
  int TIMES;
  Time_Control * Time;
  double Time_Current, value;

  Time         = Table->T;
  TIMES        = Table->T->I_Time;

  /* Each stochastic realization may be saved in a different file */
  // file[0]='\0';
  // fitxer(file, "ODE_output_", 0, ".dat"); FP = fopen(file, "w");

  /* BEGIN : Initial Conditions * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  // printf(" Before Initial_Conditions_Numerical_Integration (...)\n");
  Initial_Conditions_Numerical_Integration( Table, Table->Vector_Model_Variables );
  // printf(" After Initial_Conditions_Numerical_Integration (...). Initial Conditions:  ");

  for (k=0; k < Table->MODEL_STATE_VARIABLES; k++) {
       // printf("y[%d] = %g\t", k, Table->Vector_Model_Variables[k]);
    Table->Vector_Model_Variables_Time_0[k] = Table->Vector_Model_Variables[k];
  }
     // printf("\n");
     // Press_Key();

  Time_Current = Time->Time_Vector[0];
  if (Time->TYPE_of_TIME_DEPENDENCE > 0) Time_Dependence_Apply( Table, Time_Current );

#if defined CPGPLOT_REPRESENTATION
  Table->CPG->x_Time[0]      = Time->Time_Vector[0];
#endif

  for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
    kk = Table->OUTPUT_VARIABLE_INDEX[k];
    value = definition_OutPut_Variables(kk, Table->Vector_Model_Variables, Time->Time_Vector[0], Table);
#if defined CPGPLOT_REPRESENTATION
    Table->CPG->y_Time[k][0] = value;
#endif
    Table->Vector_Output_Variables[k]    = value;
    Table->Matrix_Output_Variables[k][0] = value;
  }
  /*   END : Initial Conditions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#if defined VERBOSE
    printf("Model Variables at time t=%g\n", Time_Current);
    char * Label = (char *)calloc( 20, sizeof(char) );
    for (k=0; k < Table->MODEL_STATE_VARIABLES; k++) {
      AssignLabel_to_Model_Variables(k, Label, Table);
      printf("y[%s] = %g\n", Label, Table->Vector_Model_Variables[k]);
    }
    printf("\n");
    free(Label);
#endif
  
#if defined VERBOSE
  printf("Output Variables at t=%g", Time_Current);
  for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
    printf("\ty[%d]=%g", k, Table->Vector_Output_Variables[k]);
  }
  printf("\n\n");
#endif


  /* BEGIN : Writing a costumized file ... */
    //fprintf(FP,"%g", Time_Current);
    // fprintf(FP,"%g", Time_Current/360.0);
    // for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
    // fprintf(FP,"\t%g", Table->Vector_Output_Variables[k]);
    // }
    // fprintf(FP,"\n");
  /*   END: Writing costumized file        */

    //printf(" Initiating Numerical Integration\n");
  for( j = 1; j < TIMES; j++ ) {
    /* This loop advances the system sequentially from
       intitial time 0 to 1st time , ...,  from time (j-1) to j, and so on.
       Note: When the system is frozen (FROZEN_SYSTEM = 1), then
       this loop does not advance the system any more
    */
    if (Table->T->TYPE_of_TIME_DEPENDENCE > 0) Time_Dependence_Apply( Table, Time_Current );
/*-------------------------------------------------------------------*/
/* B E G I N :
 *     CENTRAL POINT HERE: Numerical Integration up to the next time
 */
    State = numerical_Integration_Driver( Table, j, &Time_Current );
/*     E N D
 *-------------------------------------------------------------------*/

    if (State != GSL_SUCCESS) break;

#if defined CPGPLOT_REPRESENTATION
    Table->CPG->x_Time[j] = Time_Current;
#endif
    for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
      kk = Table->OUTPUT_VARIABLE_INDEX[k];
      value = definition_OutPut_Variables(kk, Table->Vector_Model_Variables,
					  Time_Current, Table);
#if defined CPGPLOT_REPRESENTATION
      Table->CPG->y_Time[k][j] = value;
#endif
      Table->Vector_Output_Variables[k]    = value;
      Table->Matrix_Output_Variables[k][j] = value;
    }
#if defined CPGPLOT_REPRESENTATION
      /* This should be only activated in case we want to animate ODE time evolution by
	 representing the solution as time progresses                                       */
      // C_P_G___S_U_B___P_L_O_T_T_I_N_G ( Table, j, Table->CPG->x_Time, Table->CPG->y_Time );
      // C_P_G___P_H_A_S_E____D_I_A_G_R_A_M ( Table, 0, 1, j,
      //                                      Table->CPG->y_Time );
#endif
    /* BEGIN : Writing a costumized file ... */
      // fprintf(FP,"%g", Time_Current/360.0);
      // fprintf(FP,"%g", Time_Current);
      // for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
      // fprintf(FP,"\t%g", Table->Vector_Output_Variables[k]);
      // }
      // fprintf(FP,"\n");
#if defined VERBOSE
    printf("Model Variables at time t=%g\n", Time_Current);
    char * Label = (char *)calloc( 20, sizeof(char) );
    for (k=0; k < Table->MODEL_STATE_VARIABLES; k++) {
      AssignLabel_to_Model_Variables(k, Label, Table);
      printf("y[%s] = %g\n", Label, Table->Vector_Model_Variables[k]);
    }
    printf("\n");
    free(Label);
#endif
    /*   END: Writing costumized file        */

#if defined VERBOSE
      printf("Output variables at time t = %g\n", Time_Current);
      for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
	printf("\ty[%d]=%g", k,Table->Vector_Output_Variables[k]);
      }
      printf("\n\n"); Press_Key();
#endif

   // double No_of_XX, No_of_WW, No_of_TOTAL_Women;
	 // No_of_WW = total_Female_Sexual_Workers (Table->Vector_Model_Variables, Table);
	 // No_of_XX = total_Female_Non_Sexual_Workers (Table->Vector_Model_Variables, Table);
	 // No_of_TOTAL_Women = total_Females (Table->Vector_Model_Variables, Table);
	 // printf("t=%g\tXX = %g\tWW = %g\tTotal=%g\n",Time_Current, No_of_XX, No_of_WW, No_of_TOTAL_Women);
   // Press_Key();

  }/* ------> go further to the next time step */


  // fclose(FP);

  if (State != GSL_SUCCESS) printf(" Numerical integration failed at the %dth time\n", j);

  return(State);
}
