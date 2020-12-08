#include <MODEL.h>

extern gsl_rng * r;

void Random_Initial_Condition(Parameter_Table * Table, double * Y)
{
  int i;

  double Value_0, Value_1;
  /* Definition of the state vector numerical order, from 0 to K, of model variables */
#include <Model_Variables_Code.Include.c>

  for( i = 0; i <= K; i++ ) {
    Initial_Condition_Boundary_Values(i, &Value_0, &Value_1, Table );
    Y[i] = Value_0 + gsl_rng_uniform(r) * (Value_1 - Value_0);
  }

  if(Table->RESCALING_INITIAL_TOTAL_POPULATION == 1)
    R_E_S_C_A_L_I_N_G___I_N_I_T_I_A_L___C_O_N_D_I_T_I_O_N_S ( Table, Y );
}
