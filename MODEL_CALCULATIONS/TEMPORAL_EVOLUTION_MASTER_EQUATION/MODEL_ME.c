#include "../../Include/MODEL.h"

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

  /* BEGIN : -------------------------------------------------------------------------
   * Definition Initial Condition 
   */
  // double p_1;          /* -Hp1 */ /* Resource Carrying Capacity Fraction */ 
  // double p_2;          /* -Hp2 */ /*  */ 
  Table->TOTAL_No_of_RESOURCES  = (int)(Table->p_1 * (double)Table->K_R);
  Table->TOTAL_No_of_CONSUMERS  = 200; // = Table->No_of_INDIVIDUALS; 

  assert(Table->p_2 == 1.0); 

  Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0 = (int)(Table->p_2*(double)Table->TOTAL_No_of_CONSUMERS);
  Table->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0 = Table->TOTAL_No_of_CONSUMERS - Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0;
  /* END ----------------------------------------------------------------------------
   */
  
  int No_of_CONFIGURATIONAL_STATES;
  int n_DIMENSION;
  int n_x, n_y, n_z;
  Model_Parameters_Master_Equation(Table,
				   &No_of_CONFIGURATIONAL_STATES,
				   &n_DIMENSION;
				   &n_x, &n_y, &n_z);
  
  Master_Equation * MEq = (Master_Equation *)calloc( 1, sizeof(Master_Equation) );
  Master_Equation_Allocation ( MEq,
			       No_of_CONFIGURATIONAL_STATES,
			       n_DIMENSION,
			       n_x, n_y, n_z);
  Master_Equation_Initialization ( MEq,
				   No_of_CONFIGURATIONAL_STATES,
				   n_DIMENSION,
				   n_x, n_y, n_z);
  Table->MEq = MEq; 
  
  /* Master Equation Numerical Integration                 */
  /* BEGIN: Core part (integration of the master equation) */
  printf("Entering Numerical Integration of the Master Equation...\n");   Press_Key();
  int ME_SYSTEM = master_equation_time_dynamics( Table );
  /*   END: ------------------------------------------ */

  Master_Equation_Free ( MEq );

  return(0);
}
