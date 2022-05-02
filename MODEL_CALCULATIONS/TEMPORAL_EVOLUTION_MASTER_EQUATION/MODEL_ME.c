#include <MODEL.h>

extern gsl_rng * r;

extern int TYPE_of_TIME_DEPENDENCE;

int M_O_D_E_L___M_E( Parameter_Table * Table )
{
  int i,j,k, n;
  int I_Time, no_Patch;
  int Bad_Times;
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
  
  Labels_for_Marginal_Probabilities( Table ); 
  
  /* Master Equation Numerical Integration                 */
  /* BEGIN: Core part (integration of the master equation) */
  printf("Entering Numerical Integration of the Master Equation...\n");   Press_Key();
  int ME_SYSTEM = master_equation_time_dynamics( Table );
  printf("Numerical Integration of the Master Equation succeeded!!!\n");  Press_Key();
  /*   END: ------------------------------------------ */
  
#if defined CPGPLOT_REPRESENTATION
  //  Parameter Table dependent costumized plotting is defined in
  //  ~/CPGPLOT/CPGPLOT_Parameter_Table/ files
  int TIMES           = Table->T->I_Time;
  int Input_Parameter = 0; /* The value of this model parameter appears in the title */
  //  Axes redefinition:
  Table->CPG->CPG_RANGE_X_0 = 0.0;  Table->CPG->CPG_RANGE_X_1 = Table->CPG->x_Time[TIMES-1];  
  Table->CPG->CPG_RANGE_Y_0 = 0.0;  Table->CPG->CPG_RANGE_Y_1 = (double)Table->No_of_INDIVIDUALS;
  C_P_G___S_U_B___P_L_O_T_T_I_N_G___C_U_S_T_O_M_I_Z_E_D___T_I_T_L_E(Table,
								    TIMES,
								    Table->CPG->x_Time,
								    Table->Matrix_Output_Variables,
								    Input_Parameter );
  Press_Key(); 
  C_P_G___S_U_B___P_L_O_T_T_I_N_G___C_U_S_T_O_M_I_Z_E_D___T_I_T_L_E(Table,
								    TIMES,
								    Table->CPG->x_Time,
								    Table->CPG->y_Time,
								    Input_Parameter );
#endif

  FILE * fp = fopen("Data_ODE_vs_ME_Marginals_0.dat", "w");
  for(i=0; i<TIMES; i++)
    fprintf(fp, "%g\t%g\t%g\n", Table->CPG->x_Time[i],  Table->CPG->y_Time[0][i], Table->Matrix_Output_Variables[0][i]);  
  fclose(fp);
         fp = fopen("Data_ODE_vs_ME_Marginals_1.dat", "w");
  for(i=0; i<TIMES; i++)
    fprintf(fp, "%g\t%g\t%g\n", Table->CPG->x_Time[i],  Table->CPG->y_Time[1][i], Table->Matrix_Output_Variables[1][i]);  
  fclose(fp);
  
  Master_Equation_Free ( MEq );

  return(0);
}
