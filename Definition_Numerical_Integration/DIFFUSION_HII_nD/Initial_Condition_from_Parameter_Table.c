#include <MODEL.h>

  extern gsl_rng * r;

void Initial_Condition_from_Parameter_Table(Parameter_Table * Table, double *Y)
{
  int J;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
#include <Model_Variables_Code.Include.c>

  assert(K == Table->No_of_CELLS*Table->LOCAL_STATE_VARIABLES-1);

  for (J=0; J<Table->MODEL_STATE_VARIABLES; J++)
    Y[J] = Table->Vector_Model_Variables_Time_0[J];

  if(Table->RESCALING_INITIAL_TOTAL_POPULATION == 1)
    R_E_S_C_A_L_I_N_G___I_N_I_T_I_A_L___C_O_N_D_I_T_I_O_N_S ( Table, Y );

}

void R_E_S_C_A_L_I_N_G___I_N_I_T_I_A_L___C_O_N_D_I_T_I_O_N_S ( Parameter_Table * Table, double * Y)
{
    int i;
    double g_H, E_x, x, y;
    int J;

    /* Definition of the state vector numerical order, from 0 to K, of model variables */
#include <Model_Variables_Code.Include.c>

    y = 0.0;
    for( i = 0; i <= K; i++ ) y += Y[i];
    for( i = 0; i <= K; i++ ) Y[i] = Y[i]/y * Table->INITIAL_TOTAL_POPULATION;

    /* Rescaled (Corrected) Initial Conditions */
    for (J=0; J <= K; J++)
      Y[J] = Table->Vector_Model_Variables_Time_0[J] = Y[J];
}

void Initial_Condition_Centered_into_Parameter_Table (Parameter_Table * Table, double Value)
{
  /* Initial conditions from empirical data at the initial time ( -xn 0 ) 
     Value: Number of free consumer at time 0 */
  double x; 

  int J,n,m;
  int J_X, J_Y;
  int N_X, N_Y;

  N_X = Table->No_of_CELLS_X;
  N_Y = Table->No_of_CELLS_Y;

  x = (double)Table->TOTAL_No_of_CONSUMERS; 

  if( Table->No_of_CELLS > 1 ) { 
    for (J=0; J<Table->No_of_CELLS; J++) {
      
      J_X = J/N_X;
      J_Y = J%N_X;
      
    /* In the central cell, handling consumers are equally distributed across the different handling types */
      if (N_X%2 == 0 && N_Y%2 == 0) {
	      if ( J_X == N_X/2 && J_Y == N_Y/2 ) {
          /* The central cell... */
	        for (n=0; n < Table->LOCAL_STATE_VARIABLES; n++) { 
	          m = Table->LOCAL_STATE_VARIABLES*J + n;
	          Table->Vector_Model_Variables_Time_0[m] = (x - Value)/(double)Table->LOCAL_STATE_VARIABLES;
          }
        }
        else {
          /* Out of the central cell, no handling consumers at all */
          for (n=0; n < Table->LOCAL_STATE_VARIABLES; n++) { 
	          m = Table->LOCAL_STATE_VARIABLES*J + n;
	          Table->Vector_Model_Variables_Time_0[m] = 0.0; 
          }
        }
      }
      else {
	      printf(" N_X = %d\tN_Y = %d, but both N_X and N_Y must be even!!!\n", N_X, N_Y);
	      printf(" The program will exit\n");
	      exit(0);
      }
    }
  }
  else {
    /* The system has only one cell (no spatial structure) */
    assert(Table->No_of_RESOURCES > 1 && Table->No_of_CELLS == 1);

    int Sp = (double)Table->No_of_RESOURCES;
    assert(Sp == Table->LOCAL_STATE_VARIABLES); 
   
    for (n=0; n < Table->LOCAL_STATE_VARIABLES; n++) { 
      Table->Vector_Model_Variables_Time_0[n] = (x - Value)/(double)Sp;
    }
  
    printf ("TOTAL No of CONSUMERS: %d\n", Table->TOTAL_No_of_CONSUMERS); 
    printf ("TOTAL No of FREE CONSUMERS at TIME_0: %g\n", Value);   
    assert(x -Value > 0.0);
    printf ("TOTAL No of HANDLING CONSUMERS at TIME_0: %g\n", x-Value) ; 
  /* Handling consumers are equally distributed !!! */  
  }
}

void Initial_Condition_One_Single_Cell_into_Parameter_Table (Parameter_Table * Table,
							                                               double Value_0,
							                                               double Value_1)
/* Value_0 is usually Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0,
	 Value_1 is usally Table->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0, 
   but could take other values too 
*/
{
  /* Initial conditions from empirical data at the initial time ( -xn 0 ) */

  int J,n,m;
  double x; 

  /* The system has only one cell (no spatial structure) */
  assert(Table->No_of_CELLS == 1 && Table->No_of_RESOURCES > 1);

  int Sp = (double)Table->No_of_RESOURCES;
  assert(Sp == Table->LOCAL_STATE_VARIABLES); 
  
  x = 0.0; 
  for (n=0; n < Table->LOCAL_STATE_VARIABLES; n++) { 
      Table->Vector_Model_Variables_Time_0[n] = Value_1/(double)Sp;
      x += Table->Vector_Model_Variables_Time_0[n]; 
  }
  
  printf ("TOTAL No of CONSUMERS: %d\n", Table->TOTAL_No_of_CONSUMERS); 
  printf ("TOTAL No of HANDLING CONSUMERS at TIME_0: %g\n", x); 
  x = (double)Table->TOTAL_No_of_CONSUMERS - x;  
  printf ("TOTAL No of FREE CONSUMERS at TIME_0: %g (%g)\n", x, Value_0) ; 
  // printf ("Handling consumers are equally distributed accross behavioural types!!!\n");
}

void Initial_Condition_All_Patches_the_Same_into_Parameter_Table (Parameter_Table * Table, 
                                                                  double Value)
{
  /* Initial conditions from empirical data at the initial time ( -xn 0 ) */

  int J,n,m; 
  int N_X, N_Y;

  N_X = Table->No_of_CELLS_X;
  N_Y = Table->No_of_CELLS_Y;

  for (J=0; J<Table->No_of_CELLS; J++) {
    for (n=0; n < Table->LOCAL_STATE_VARIABLES; n++) {
      m = Table->LOCAL_STATE_VARIABLES*J + n;
      Table->Vector_Model_Variables_Time_0[m] = Value;
    }
  }
}

  void Common_Initial_Condition_Command_Line_Arguments_into_Table(Parameter_Table *Table)
{
  /* BEGIN : -------------------------------------------------------------------------
   * Definition Initial Condition:  
   */
  /* This definition is contingent to TYPE of MODEL at work from the pre-defined family 
	 * of models:
	 * DIFFUSION_BD_2D, DIFFUSION_HII_1D, DIFFUSION_BD_2D, and 
   *  DIFFUSION_HII_nD (so far)
	 * ( Models where the TOTAL_No_of_CONSUMERS is a CONSTANT )
   */
  /* double p_1;         */ /* -Hp1 */ /* Resource Carrying Capacity Fraction */ 
  /* double p_2;         */ /* -Hp2 */ /* See below the definition of the     */
                                      /* TOTAL_No_of_FREE_CONSUMERS_TIME_0   */
  Table->TOTAL_No_of_RESOURCES  = (int)(Table->p_1 * (double)Table->K_R);
  Table->TOTAL_No_of_CONSUMERS  = Table->No_of_INDIVIDUALS;  /* -HN 20 as input argument */ 

  assert(Table->p_2 <= 1.0 && Table->p_2 >= 0.0);  // 
  assert(Table->p_1 <= 1.0 && Table->p_1 >= 0.0);  // Fractions!!!  
  
  Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0     = (int)(Table->p_2*(double)Table->TOTAL_No_of_CONSUMERS);
  Table->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0 = Table->TOTAL_No_of_CONSUMERS - Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0;
  Table->TOTAL_No_of_FREE_CONSUMERS            = (int)(Table->p_2*(double)Table->TOTAL_No_of_CONSUMERS);
  Table->TOTAL_No_of_HANDLING_CONSUMERS        = Table->TOTAL_No_of_CONSUMERS - Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0;
  /* END ----------------------------------------------------------------------------
	 */
}