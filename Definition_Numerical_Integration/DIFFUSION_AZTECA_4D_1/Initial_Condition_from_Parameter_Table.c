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

void Initial_Condition_Centered_into_Parameter_Table (Parameter_Table * Table,
						                                          double Value)
{
  /* Initial conditions from empirical data at the initial time ( -xn 0 ) */
  int No_of_CHOICES;
  int J_X_R, J_Y_R;
  int J,n,m,k;
  int J_X, J_Y;
  int N_X, N_Y;

  N_X = Table->No_of_CELLS_X;
  N_Y = Table->No_of_CELLS_Y;

  if( Table->No_of_CELLS > 1 ) { 
    for (J=0; J<Table->No_of_CELLS; J++) {
      J_X = J/N_X;
      J_Y = J%N_X;
      
      /* Some Workers (0) at the center */
      if (N_X%2 == 0 && N_Y%2 == 0) {
	      if ( J_X == N_X/2 && J_Y == N_Y/2 ) {
	        n=0;                                             /* n=0: Workers */
	        m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = Value;
          n=1;                                             /* n=1: Queens */
	        m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = Value;
          n=2;                                             /* n=2: Flies */
	        m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = 0.0;
          n=3;                                             /* n=3: Parasitized Workers */
	        m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = 0.0;   
      	}
	      else {
          /* Workers (F, 0) */
          n = 0; 
	        m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = 0.0;
          /* Queens (Q, 0) */
          n = 1;
          m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = 0.0;
          /* Flies (F, 2) */
          n = 2;
          m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = 0.0;
          /* Larvae](WF, 3) (0) */
          n = 3; 
          m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = 0.0;
	      }
      }
      else {
	      printf(" N_X = %d\tN_Y = %d, but both N_X and N_Y must be even!!!\n",
	             N_X, N_Y);
	      printf(" The program will exit\n");
	      exit(0);
      }
    }
    
    /* Number of Cells to be chosen randomly */
    // No_of_CHOICES = (int)(Table->p_1 * Table->No_of_CELLS);
    No_of_CHOICES = 3; 
    assert(Table->No_of_CELLS > No_of_CHOICES);


    for (k = 0; k<No_of_CHOICES; k++){    
      /* All other Species at scattered at a number of randomly chosen cells 
         located at random placements (J_X_R, J_Y_R) */
      J_X_R = (int)(Table->No_of_CELLS_X * gsl_rng_uniform(r));
      J_Y_R = (int)(Table->No_of_CELLS_Y * gsl_rng_uniform(r));
      
        /* Flies (2), Workers (0), and Queens (1) scattered at random 
           in No_of_CHOICES cells 
        */  
      for (J=0; J<Table->No_of_CELLS; J++) {
	      J_X = J/N_X;
	      J_Y = J%N_X;
	
	      if ( J_X == J_X_R && J_Y == J_Y_R ) {
          /* Workers (F, 0) */
          n = 0; 
	        m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = Value;
          /* Queens (Q, 0) */
          n = 1;
          m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = Value;
          /* Flies (F, 2) */
          n = 2;
          m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = Value;
          /* Larvae](WF, 3) (0) */
          n = 3; 
          m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = 0.0;
	      }
	      else {
          /* Workers (F, 0) */
          n = 0; 
	        m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = 0.0;
          /* Queens (Q, 0) */
          n = 1;
          m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = 0.0;
          /* Flies (F, 2) */
          n = 2;
          m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = 0.0;
          /* Larvae](WF, 3) (0) */
          n = 3; 
          m = Table->LOCAL_STATE_VARIABLES*J + n;
	        Table->Vector_Model_Variables_Time_0[m] = 0.0;
	      }
      }
    }
  }
  else {
    /* The system has only one cell (no spatial structure) */
    assert(Table->No_of_RESOURCES == 1 && Table->No_of_CELLS == 1);

    Table->Vector_Model_Variables_Time_0[0] = Value;    /* Initial value for workers         */
    Table->Vector_Model_Variables_Time_0[1] = Value;  /* Initial value for queens */
    Table->Vector_Model_Variables_Time_0[2] = Value; /* Initial value for flies */
    Table->Vector_Model_Variables_Time_0[3] = 0.0;   /* Initial value for parasitized workers */
  }
}

void Initial_Condition_One_Single_Cell_into_Parameter_Table (Parameter_Table * Table,
							                                               double Value_0,
							                                               double Value_1)
{
  /* Initial conditions from empirical data at the initial time ( -xn 0 ) */

  int J,n,m;
 
  /* The system has only one cell (no spatial structure) */
  assert(Table->No_of_CELLS     == 1);
  assert(Table->No_of_RESOURCES == 1);

  Table->Vector_Model_Variables_Time_0[0] = Value_0; /* Initial value for workers            */
  Table->Vector_Model_Variables_Time_0[1] = Value_0; /* Initial value for queens             */
  Table->Vector_Model_Variables_Time_0[2] = Value_1; /* Initial value for free flies         */
  Table->Vector_Model_Variables_Time_0[3] = 0.0;     /* Initial value for fly larvae (WF)    */
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
    /* Resources */
    n = 0;  
    m = Table->LOCAL_STATE_VARIABLES*J + n;
    Table->Vector_Model_Variables_Time_0[m] = Value;

    /* Consumers */
    n = 1;  
    m = Table->LOCAL_STATE_VARIABLES*J + n;
    Table->Vector_Model_Variables_Time_0[m] = Value;
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
	 * This initial Condition involves no triplets at time t = 0.0 because the sum of states
	 * should add up the TOTAL No of CONSUMERS 
   */
}
