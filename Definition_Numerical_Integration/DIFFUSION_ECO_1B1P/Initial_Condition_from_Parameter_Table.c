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
  int Sp; 

  N_X = Table->No_of_CELLS_X;
  N_Y = Table->No_of_CELLS_Y;

  if( Table->No_of_CELLS > 1 ) { 
    for (J=0; J<Table->No_of_CELLS; J++) {
      J_X = J/N_X;
      J_Y = J%N_X;
      
      /* Species 0 at the center */
      /* Species 1 (or any other bacterial type) at center */       
      if (N_X%2 == 0 && N_Y%2 == 0) {
	      if ( J_X == N_X/2 && J_Y == N_Y/2 ) {
	        Sp=0;                                            
	        m = Table->LOCAL_STATE_VARIABLES*J + Sp;    
	        Table->Vector_Model_Variables_Time_0[m] = Value;

          Sp=1;                                            
	        m = Table->LOCAL_STATE_VARIABLES*J + Sp;    
	        Table->Vector_Model_Variables_Time_0[m] = Value;
      	}
	      else {
          Sp=0;                                            
	        m = Table->LOCAL_STATE_VARIABLES*J + Sp;    
	        Table->Vector_Model_Variables_Time_0[m] = 0.0;

          Sp=1;                                            
	        m = Table->LOCAL_STATE_VARIABLES*J + Sp;    
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

    No_of_CHOICES = 4; /* Number of Cells to be chosen */
    for (k = 0; k<No_of_CHOICES; k++){    
      /* All other Species at scattered at a number of randomly chosen cells 
         located at random placements (J_X_R, J_Y_R) */
      J_X_R = (int)(Table->No_of_CELLS_X * gsl_rng_uniform(r));
      J_Y_R = (int)(Table->No_of_CELLS_Y * gsl_rng_uniform(r));
      
      n = (int)(Table->No_of_RESOURCES * gsl_rng_uniform(r));   /* Ten choices of a species type 
                                                                   at random (in B1P1, this is 0 or 1) 
                                                                */
      assert(Table->No_of_RESOURCES ==  Table->LOCAL_STATE_VARIABLES);

      for (J=0; J<Table->No_of_CELLS; J++) {
        J_X = J/N_X;
        J_Y = J%N_X;

        if ( J_X == J_X_R && J_Y == J_Y_R ) {
  
          m = Table->LOCAL_STATE_VARIABLES*J + n;
          Table->Vector_Model_Variables_Time_0[m] = Value;
        }
        else {
  
          m = Table->LOCAL_STATE_VARIABLES*J + n;
          Table->Vector_Model_Variables_Time_0[m] = 0.0;
        }
      }
    }
  }
  else {
    /* The system has only one cell (no spatial structure) */
    assert(Table->No_of_RESOURCES == 2 && Table->No_of_CELLS == 1);

    Sp = (int)(Table->No_of_RESOURCES * gsl_rng_uniform(r)); 
    Table->Vector_Model_Variables_Time_0[Sp] = Value; /* Initial value for strain Sp */

    assert(Table->No_of_RESOURCES > Sp);

    /* The rest of values are zero because calloc definition of Vector_Model_Variables_Time_0                */
  }
}

void Initial_Condition_One_Single_Cell_into_Parameter_Table (Parameter_Table * Table,
							                                               double Value_0,
							                                               double Value_1)
{
  /* Initial conditions from empirical data at the initial time ( -xn 0 ) */
  int j; 
  int R, Sp_0, Sp_1;

  /* The system has only one cell (no spatial structure) */
  assert(Table->No_of_CELLS == 1);

  /* The choice of two random types will be initiated with non-zero populations: */
  Sp_0 = 0;
  Sp_1 = 1;
   
  for(j=0; j<Table->No_of_CELLS; j++) {
           
      R  = j * Table->LOCAL_STATE_VARIABLES + Sp_0;
      Table->Vector_Model_Variables_Time_0[R]  = Value_0; /* Initial value for resources Species Sp_0  */

      R  = j * Table->LOCAL_STATE_VARIABLES + Sp_1;
      Table->Vector_Model_Variables_Time_0[R]  = Value_1; /* Initial value for resources Species Sp_1  */
        
  }  
}    


void Initial_Condition_All_Patches_the_Same_into_Parameter_Table (Parameter_Table * Table,
								                                                  double Value)
{
  /* Initial conditions from empirical data at the initial time ( -xn 0 ). All types take
     the same intial value across cells in the system. */

  int J,m, Sp;  

  for (J=0; J<Table->No_of_CELLS; J++) {
    for(Sp = 0; Sp < Table->No_of_RESOURCES; Sp++) {
    /* Resources Propagules Sp */
      m = Table->LOCAL_STATE_VARIABLES*J + Sp ;
      Table->Vector_Model_Variables_Time_0[m] = Value;
    }
  }
}

void Common_Initial_Condition_Command_Line_Arguments_into_Table(Parameter_Table *Table)
{
  /* BEGIN : -------------------------------------------------------------------------
   * Definition Initial Condition:  
   */
  
  assert(Table->Chi_C_0 <= 1.0 && Table->Chi_C_0 >= 0.0);    /* Conjugation probability (Coded as Xhi_C_0, but labeled as Xhi) */   
  
  assert(Table->p_2 <= 1.0 && Table->p_2 >= 0.0);            /* Competition-induced Death Probability  (Coded as p_2, but labeled as p_C) */    
  
  assert(Table->p_1 <= 1.0 && Table->p_1 >= 0.0);            /* Segregation error at reproduction (Coded as p_1, but labeled as \epsilon) */    
  
 /* END ----------------------------------------------------------------------------
  */
}
