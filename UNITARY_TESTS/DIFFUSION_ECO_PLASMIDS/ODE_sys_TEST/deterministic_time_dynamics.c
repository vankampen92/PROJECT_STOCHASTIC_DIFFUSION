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
#define BIO_NEGATIVE_VALUE -1.0E-10

int Deterministic_Time_Dynamics( Parameter_Table * Table )
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
  */
  int NEGATIVE_VALUE; 
  int i; int State;
  
  int j, k, kk;
  int TIMES;
  Time_Control * Time;
  double Time_Current, value;

  Time         = Table->T;
  TIMES        = Table->T->I_Time;

  /* BEGIN : Initial Conditions * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  printf(" Before Initial_Conditions_Numerical_Integration (...)\n");
  Initial_Conditions_Numerical_Integration( Table, Table->Vector_Model_Variables );
  printf(" After Initial_Conditions_Numerical_Integration (...). Initial Conditions:  ");

  assert( Table->LOCAL_STATE_VARIABLES == Table->No_of_RESOURCES );
  assert( Table->MODEL_STATE_VARIABLES == Table->No_of_RESOURCES );  
  assert( Table->No_of_CELLS == 1);

  Time_Current = Time->Time_Vector[0];

#if defined CPGPLOT_REPRESENTATION
    Table->CPG->x_Time[0]      = Time->Time_Vector[0];
    for(k=0; k < Table->No_of_RESOURCES; k++){
      Table->CPG->y_Time[k][0] = Table->Vector_Model_Variables[k];
    }
#endif
  /*   END : Initial Conditions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  printf(" Initiating Numerical Integration\n");
  for( j = 1; j < TIMES; j++ ) {
    /* This loop advances the system sequentially from
       intitial time 0 to 1st time , ...,  from time (j-1) to j, and so on.
       Note: When the system is frozen (FROZEN_SYSTEM = 1), then
       this loop does not advance the system any more
    */
    /*-------------------------------------------------------------------*/
    /* B E G I N :
     *     CENTRAL POINT HERE: Numerical Integration up to the next time
     */
    State = numerical_Integration_Driver( Table, j, &Time_Current );
    /*     E N D
     * ------------------------------------------------------------------*/

    if (State != GSL_SUCCESS) break;

#if defined CPGPLOT_REPRESENTATION
    Table->CPG->x_Time[j] = Time_Current;

    for(k=0; k < Table->No_of_RESOURCES; k++){
      Table->CPG->y_Time[k][j] = Table->Vector_Model_Variables[k];
    }
#endif
    
#if defined CPGPLOT_REPRESENTATION
    /* This should be only activated in case we want to animate ODE time evolution by
	     representing the solution as time progresses                                       */
    Table->CPG->CPG_RANGE_X_0 = Table->T->Time_0; 
    Table->CPG->CPG_RANGE_X_1 = Table->T->Time_1;
    
    Table->CPG->CPG_RANGE_Y_0 = 0.0;              
    Table->CPG->CPG_RANGE_Y_1 = (double)Table->K_R;

    int TimeEvoPlot =  CPGPLOT___X_Y_n___P_L_O_T_T_I_N_G( Table->CPG,
                                                          j, Table->No_of_RESOURCES,
                                                          Table->CPG->x_Time, Table->CPG->y_Time,
                                                          "Time", "Abundance", "",
                                                          1, 1 );             
    getchar();
#endif
  }/* ------> go further to the next time step */

  return(State);
}
