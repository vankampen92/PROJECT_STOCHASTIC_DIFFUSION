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
     avancing from a time to the next as stored in Time->Time_Vector[], and save a file 
     corresponding to this numerical integration.

     This function makes two essential calls:

     1. Initial_Condition_Master_Equation (...), which sets up initial conditions. 

     2. master_equation_driver(...), which performs the actual numerical integration

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

  // file[0]='\0';
  // fitxer(file, "ME_output_", 0, ".dat"); FP = fopen(file, "w");

  /* BEGIN : Initial Conditions * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#if defined VERBOSE
  printf(" Before Initial Condition Master Equation (...)\n");
#endif
  
  Initial_Condition_Master_Equation( Table, Table->MEq->Probability_Distribution );

#if defined VERBOSE
  printf(" After Initial Condition Master Equation (...).");
#endif
  
  for (k=0; k < Table->MEq->No_of_CONFIGURATIONAL_STATES; k++) {
    Table->MEq->Probability_Distribution_Time_0[k] = Table->MEq->Probability_Distribution[k];
  }
  
  Time_Current = Time->Time_Vector[0];
  if (Time->TYPE_of_TIME_DEPENDENCE > 0) Time_Dependence_Apply( Table, Time_Current );
  /*   END : Initial Conditions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#if defined VERBOSE
  printf(" Initiating Numerical Integration\n");
#endif
  for( j = 1; j < TIMES; j++ ) {
    /* This loop advances the system sequentially from
       intitial time 0 to 1st time , ...,  from time (j-1) to j, and so on.
       Note: When the system is frozen (FROZEN_SYSTEM = 1), then
       this loop does not advance the system any more
    */
    if (Table->T->TYPE_of_TIME_DEPENDENCE > 0) Time_Dependence_Apply( Table, Time_Current );
/*-------------------------------------------------------------------*/
/* B E G I N :
 *  CENTRAL POINT HERE: Numerical Integration of the master equation up to the next time
 */
    State = master_equation_driver( Table, j, &Time_Current );
    /* Normalization: */
    Normalization_Master_Equation(Table->MEq->Probability_Distribution,
				  Table->MEq->No_of_CONFIGURATIONAL_STATES);
    
/*     E N D : ------------------------------------------------------*/

    if (State != GSL_SUCCESS) break;

#if defined CPGPLOT_REPRESENTATION
    /* This should be only activated in case we want to animate ODE time evolution by
       representing the solution as time progresses                                       
    */
    
    C_P_G___P_R_O_B_A_B_I_L_I_T_Y___D_I_S_T_R_I_B_U_T_I_O_N ( Table, j, Time_Current );
#endif

  }/* ------> go further to the next time step */

#if defined VERBOSE
  printf(" Numerical Integration of the Master Equaiton ended!!!\n");
#endif

  // fclose(FP);

#if defined VERBOSE
  if (State != GSL_SUCCESS ) {
    printf(" Numerical integration failed at the %dth time\n", j);
  }
#endif
  
  return(State);
}

