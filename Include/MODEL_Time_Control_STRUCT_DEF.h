typedef struct Time_Controlinfo
{
  #include "include.Time_Control.global.h" /* Control variable I_Time */

  double * Time_Vector;

  double ** summ;
  double ** summ_var;
  int * count;

  double ** AVE;
  double ** VAR;
  double * time_DEF;

  double ** Accumulated_Variable;

  Stochastic_Rate * Rate;

  double * Vector_of_Rates; 

  double *** Variable; 

  double ** Time_Vector_Real; 
}Time_Control;

/* void T_I_M_E___S_T_R_U_C_T___A_L_L_O_C( Time_Control * Time, Parameter_Table * P, int I_Time); */

/* void  T_I_M_E___S_T_R_U_C_T___U_P_L_O_A_D( Time_Control * Time, Parameter_Table * Table, int I_Time); */

/* void T_I_M_E___S_T_R_U_C_T___F_R_E_E( Time_Control * Time, Parameter_Table * P ); */
