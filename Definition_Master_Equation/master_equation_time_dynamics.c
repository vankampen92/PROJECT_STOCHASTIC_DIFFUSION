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
  printf(" Initiating Numerical Integration Master Equaiton \n");
#endif
  Marginal_Probability_Calculation ( Table );           /* At time zero */
  Marginal_Probability_Averages_Calculation ( Table );  /* At time zero */
  Print_Probability_Distribution ( Table );
  printf("t = %g\t<n> = %g\t<m> = %g\n",
	 Time_Current,
	 Table->MEq->Vector_Model_Variables[0],
	 Table->MEq->Vector_Model_Variables[1]);  
  
#if defined CPGPLOT_REPRESENTATION
  int SAME = 0;
  assert(Table->MEq->n_DIMENSION <= 2);  
  C_P_G___M_A_R_G_I_N_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Table, j,
						      0, Time_Current,
						      SAME );
    
  C_P_G___M_A_R_G_I_N_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Table, j,
						      1, Time_Current,
						      SAME );
  Press_Key();
#endif
  // SAME = 1
  for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
    kk = Table->OUTPUT_VARIABLE_INDEX[k];
    value = definition_OutPut_Variables(kk,
					Table->MEq->Vector_Model_Variables, Time->Time_Vector[0],
					Table);
    Table->Matrix_Output_Variables[k][0] = value;
  }
  
  
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
    
/*     E N D : ------------------------------------------------------*/

    if (State != GSL_SUCCESS) break;

#if defined CPGPLOT_REPRESENTATION
    /* This should be only activated in case we want to animate ODE time evolution by
       representing the solution as time progresses                                       
    */
    Marginal_Probability_Calculation ( Table );
    Marginal_Probability_Averages_Calculation ( Table );  /* At time zero */
    for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
      kk = Table->OUTPUT_VARIABLE_INDEX[k];
      value = definition_OutPut_Variables(kk,
					  Table->MEq->Vector_Model_Variables, Time_Current, 
					  Table);
      Table->Matrix_Output_Variables[k][j] = value;
    }
    
    assert(Table->MEq->n_DIMENSION <= 2);  
    C_P_G___M_A_R_G_I_N_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Table, j,
							0, Time_Current,
							SAME );
    
    C_P_G___M_A_R_G_I_N_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Table, j,
							1, Time_Current,
							SAME );
    // SAME = 1;
    Press_Key(); 
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

void C_P_G___M_A_R_G_I_N_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Parameter_Table * Table,
							 int j,
							 int n, 
							 double Time_Current,
							 int SAME )
{
  /* This function plots marginal probability distributions. As input, it takes Table. 
     A member of Table is MEq, which stores all the necessary information to plot the 
     whole distribution and, of course, its marginals.
     
     This function will plot the marginals, depending on an input parameter, an integer,
     n, defining the marginal to plot: 
     . n=0, first dimension
     . n=1, second dimension
     . n=2, third dimension
  */
  
  /* Notice that the probability distribution and marginals in the two dimensions are 
     associated to a Time_Current around the j-th time in Time_Vector[j] 
  */
  int i;
  int No_of_POINTS;
  static double Current_Time = 0.0;
  double Last_Time;

  Parameter_CPGPLOT * CPG = Table->CPG; 
  Master_Equation * ME = Table->MEq;
  
  double * y;
  double * x; 
  
  if(n == 0)       {
    y = ME->P_n_Marginal;
    x = (double *)calloc(ME->n_x, sizeof(double));
    for(i=0; i<ME->n_x; i++) x[i] = (double)i;
    No_of_POINTS = ME->n_x;
  }
  else if (n == 1) {
    y = ME->P_m_Marginal;
    x = (double *)calloc(ME->n_y, sizeof(double));
    for(i=0; i<ME->n_y; i++) x[i] = (double)i;
    No_of_POINTS = ME->n_y;
  }
  else {
    printf("Marginal probabilities for higher dimensions (n > %d) are not coded\n", n);
    printf("The program will exit\n");
    Press_Key();
    exit(0); 
  }
  
  CPG->CPG_RANGE_X_0 = -0.5;  CPG->CPG_RANGE_X_1 = No_of_POINTS - 0.5;
  CPG->CPG_RANGE_Y_0 =  0.0;  CPG->CPG_RANGE_Y_1 = 1.0;  
  
  double x_Time_Position = 0.90 * (float)CPG->CPG_RANGE_X_1;
  double y_Time_Position = 0.90 * (float)CPG->CPG_RANGE_Y_1;

  char * Plot_Title = (char *)calloc( 100, sizeof(char));
  char * Plot_Time  = (char *)calloc( 50, sizeof(char));
  char * Time_Eraser = (char *)calloc(50, sizeof(char));
  
  Last_Time     = Current_Time;
  Current_Time  = Time_Current; 
  sprintf(Plot_Time, "Time = %5.2f", Current_Time);
  sprintf(Time_Eraser, "Time = %5.2f", Last_Time);

  Plot_Title[0] ='\0';
  sprintf(Plot_Title, "Time = %5.2f", Current_Time);
  
  CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T ( CPG, SAME,
							No_of_POINTS, 
							x, y,
							ME->Marginal_Probability_Label[n], 
							"Probability", 
							Plot_Title,
							1, 1 );
  cpgsch(2.0);
  cpgsci(0);
  cpgptxt(x_Time_Position, y_Time_Position, 0.0, 1.0, Time_Eraser);
  cpgsci(1);
  cpgptxt(x_Time_Position, y_Time_Position, 0.0, 1.0, Plot_Time);
  // cpgsch(char_Size);

  free(Plot_Title); free(Plot_Time); free(Time_Eraser);
  free(x); 
}
