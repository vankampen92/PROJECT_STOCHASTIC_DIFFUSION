#include <MODEL.h>

extern gsl_rng * r;

extern int TYPE_of_TIME_DEPENDENCE;

int M_O_D_E_L___M_E( Parameter_Table * Table )
{
  int i,j,k, n;
  int I_Time;
  double t; 
  Time_Control * Time;
  Time_Dependence_Control * TDC; 

  Time = Table->T;
  TDC  = Table->TDC; 

  int No_of_CONFIGURATIONAL_STATES;
  int n_DIMENSION;
  int n_x, n_y, n_z;

  Model_Parameters_Master_Equation(Table,
				                           &No_of_CONFIGURATIONAL_STATES,
				                           &n_DIMENSION,
				                           &n_x, &n_y, &n_z);
  
  Master_Equation * MEq = (Master_Equation *)calloc( 1, sizeof(Master_Equation) );
  /* MEq and Table will be two data structures that point to each other */
  Table->MEq = MEq;           
  MEq->Table = Table;
  
  Master_Equation_Allocation ( MEq,
			                         No_of_CONFIGURATIONAL_STATES,
			                         n_DIMENSION,
			                         n_x, n_y, n_z );

  Master_Equation_Initialization ( MEq,
				                           No_of_CONFIGURATIONAL_STATES,
				                           n_DIMENSION,
				                           n_x, n_y, n_z );                             
  
  Master_Equation_Configurational_State_Setup ( MEq );
  
  Labels_for_Marginal_Probabilities( Table ); 

  /* BEGIN : Stationary Probability Distribution (only available for some of the models) -*/
  #ifdef DIFFUSION_BD_2D
    Stationary_Probability_Distribution( Table );
    printf(" Theoretical Stationary Probability Distribution has been calculated.\n");
  #elif defined DIFFUSION_HII_1D  
    Stationary_Probability_Distribution( Table );
    printf(" Theoretical Stationary Probability Distribution has been calculated.\n");
  #elif defined DIFFUSION_HII_nD
    Stationary_Probability_Distribution( Table );
    printf(" Theoretical Stationary Probability Distribution has been calculated.\n");
    /* Only exact when all Nu's are the same!!! */
    assert_HLL_nD_Equal_Nus( Table );
    Evolving_Probability_Distribution( Table );
    printf(" Theoretical Time Evolving Probability Distribution has been calculated.\n");
  #endif
  Print_Press_Key(0, 1, ".");
  /*  END : ------------------------------------------------------------------------------*/

  /***********************************************************************************************/
  /* Master Equation Numerical Integration                 */
  /* BEGIN: Core part (integration of the master equation) */
  printf("Entering Numerical Integration of the Master Equation...\n");   Print_Press_Key(0, 1, ".");
  int ME_SYSTEM = master_equation_time_dynamics( Table );
  printf("Numerical Integration of the Master Equation succeeded!!!\n");  Print_Press_Key(0, 1, ".");
  /*   END: ------------------------------------------ */
  /***********************************************************************************************/
  
#if defined STOCHASTIC_REALIZATIONS  
 /* BEGIN : Representing distributions associated to output variables across 
            stochastic realizations 
 */  
#if defined CPGPLOT_REPRESENTATION
  int SAME_PLOT = 1;
  // j = Table->T->I_Time - 1;
  
  for ( i=0; i < Table->MEq->n_DIMENSION; i++ ) {

  #ifdef DIFFUSION_BD_2D
    C_P_G___S_T_A_T_I_O_N_A_R_Y___D_I_S_T_R_I_B_U_T_I_O_N ( Table, i,  SAME_PLOT );
  #elif defined DIFFUSION_HII_1D  
    C_P_G___S_T_A_T_I_O_N_A_R_Y___D_I_S_T_R_I_B_U_T_I_O_N ( Table, i,  SAME_PLOT );
  #elif defined DIFFUSION_HII_nD  
    C_P_G___S_T_A_T_I_O_N_A_R_Y___D_I_S_T_R_I_B_U_T_I_O_N ( Table, i,  SAME_PLOT );
  #endif
    Print_Press_Key(1,0,".");
    for(j=1; j<Table->T->I_Time; j++) {
      C_P_G___E_M_P_I_R_I_C_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Table, j, i,
  							                                            Table->T->Time_Vector[j],
  							                                            SAME_PLOT );
      Print_Press_Key(1, 1,".");
    }
  }
#endif
 /*   END : Empirical Distribuiton Representation ----------------------------------*/
  
  /* Saving also the last time in a file */
  j = Table->T->I_Time - 1;
  if(Table->MEq->n_DIMENSION <= 2)
    Saving_Empirical_Distribution_vs_ME_Numerical_Integration ( Table,
							                                                  j, 
                                                                Table->T->Time_Vector[j] );  
#endif
  
#if defined CPGPLOT_REPRESENTATION
  //  Parameter Table dependent costumized plotting is defined in
  //  ~/CPGPLOT/CPGPLOT_Parameter_Table/ files
  // Plotting the temporal evolution of the averages of the dynamical variables (which
  // should match the output variables for this work) calculated with the marginals from
  // the numerical integration of the master equation vs the densities calculated with
  // the ODE system.
  
  int TIMES           = Table->T->I_Time;
  //  Axes redefinition:
  Table->CPG->CPG_RANGE_X_0 = 0.0;  Table->CPG->CPG_RANGE_X_1 = Table->CPG->x_Time[TIMES-1];  
  Table->CPG->CPG_RANGE_Y_0 = 0.0;  Table->CPG->CPG_RANGE_Y_1 = (double)Table->No_of_INDIVIDUALS;
  Table->CPG->CPG_SCALE_Y   = 1;    Table->CPG->CPG_SCALE_X   = 1;
  SAME_PLOT = 0;
  C_P_G___S_U_B___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T( SAME_PLOT,
						       Table, TIMES,
						       Table->CPG->x_Time,
						       Table->Matrix_Output_Variables);
  SAME_PLOT = 1;
  Print_Press_Key(1,0,".");
  /* New colors and lines */
  Table->CPG->color_Index   = 3;
  Table->CPG->type_of_Width = 2;
  Table->CPG->type_of_Line  = 2;
  Table->CPG->type_of_Symbol = 5; 
  C_P_G___S_U_B___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T( SAME_PLOT,
						                                           Table, TIMES,
						                                           Table->CPG->x_Time,
						                                           Table->CPG->y_Time );
#endif

  FILE * fp = fopen("Data_ODE_vs_ME_Marginals_0.dat", "w");
  for(i=0; i<TIMES; i++)
    fprintf(fp, "%g\t%g\t%g\n", Table->CPG->x_Time[i],
	    Table->CPG->y_Time[0][i], Table->Matrix_Output_Variables[0][i]);  
  fclose(fp);
  
  if( Table->SUB_OUTPUT_VARIABLES > 1 ){
    fp = fopen("Data_ODE_vs_ME_Marginals_1.dat", "w");
    for(i=0; i<TIMES; i++)
      fprintf(fp, "%g\t%g\t%g\n", Table->CPG->x_Time[i],
	      Table->CPG->y_Time[1][i], Table->Matrix_Output_Variables[1][i]);  
    fclose(fp);
  }
  
  Master_Equation_Free ( MEq );

  return(0);
}





