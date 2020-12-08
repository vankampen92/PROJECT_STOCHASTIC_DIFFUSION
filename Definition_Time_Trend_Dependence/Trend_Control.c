#include <MODEL.h>
#include <include.Trend_Control.extern.h>

void T_R_E_N_D___C_O_N_T_R_O_L___F_R_E_E( Trend_Control * Tr )
{
  free(Tr);
}

void T_R_E_N_D___C_O_N_T_R_O_L___A_L_L_O_C( Trend_Control * Tr, Parameter_Table * P )
{
  /* */
}

void  T_R_E_N_D___C_O_N_T_R_O_L___U_P_L_O_A_D( Trend_Control * T, Parameter_Table * Table )
{

  /* Trend parameter definition */
  T->Tr_No_of_Jumps = Tr_No_of_Jumps;

  T->Tr_value_0     = Tr_value_0;
  T->Tr_value_1     = Tr_value_1;
  T->Tr_value_i     = Tr_value_i;

  T->Tr_Time_0      = Tr_Time_0;
  T->Tr_Time_1      = Tr_Time_1;
  T->Tr_Time_i      = Tr_Time_i;

  /* Default: */
  T->Tr_Input_Parameter = Tr_Input_Parameter;
  T->TYPE_of_TREND      = TYPE_of_TREND;
  /* -------- */

  Table->Tr           = T;
}

void S_E_T_U_P___T_R_E_N_D___T_R_I_A_N_G_L_E ( Trend_Control * T,
					       Time_Control * Time,
					       int Input_Parameter )
{
  double Value_0, Value_1;

  T->Tr_Time_0      = Time->Time_0;
  T->Tr_Time_1      = Time->Time_1;
  /* Note that Time_i will be here an intermediate time! */
  T->Tr_Time_i      = (Time->Time_1 - Time->Time_0)/2.0;

  Parameter_Space_Boundary_for_Trends(Input_Parameter, &Value_0, &Value_1);

  T->Tr_value_0     = Value_0;
  T->Tr_value_1     = Value_0;
  /* Note that at the intermediate time we have the highest value */
  T->Tr_value_i     = Value_1;

  T->Tr_Input_Parameter = Input_Parameter;
  T->TYPE_of_TREND         = 0;
}

void S_E_T_U_P___T_R_E_N_D___L_I_N_E_A_R ( Trend_Control * T,
					   Time_Control * Time,
					   int Input_Parameter,
					   int Signature )
{
  double Value_0, Value_1;

  T->Tr_Time_0      = Time->Time_0;
  T->Tr_Time_1      = Time->Time_1;

  Parameter_Space_Boundary_for_Trends(Input_Parameter, &Value_0, &Value_1);

  if( Signature == +1 ){
    /* Increasing trend */
    T->Tr_value_0     = Value_0;
    T->Tr_value_1     = Value_1;
  }
  else {
    /* Decreasing trend */
    T->Tr_value_0     = Value_1;
    T->Tr_value_1     = Value_0;
  }

  T->Tr_Input_Parameter = Input_Parameter;
  T->TYPE_of_TREND         = 1;
}
