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
  T->Tr_value_2     = Tr_value_2;
  T->Tr_value_3     = Tr_value_3;
  T->Tr_value_4     = Tr_value_4;
  T->Tr_value_5     = Tr_value_5;
  T->Tr_value_6     = Tr_value_6;
  T->Tr_value_7     = Tr_value_7;
  T->Tr_value_8     = Tr_value_8;
  T->Tr_value_9     = Tr_value_9;
  
  T->Tr_value_i     = Tr_value_i;

  T->Tr_Time_0      = Tr_Time_0;
  T->Tr_Time_1      = Tr_Time_1;
  T->Tr_Time_2      = Tr_Time_2;
  T->Tr_Time_3      = Tr_Time_3;
  T->Tr_Time_4      = Tr_Time_4;
  T->Tr_Time_5      = Tr_Time_5;
  T->Tr_Time_6      = Tr_Time_6;
  T->Tr_Time_7      = Tr_Time_7;
  T->Tr_Time_8      = Tr_Time_8;
  T->Tr_Time_9      = Tr_Time_9;
  
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
  T->TYPE_of_TREND         = 2;
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

  if      (Signature == +1) T->TYPE_of_TREND         = 3;
  else if (Signature == -1) T->TYPE_of_TREND         = 5;
  else {
    printf("Linear trend ill-defined. Error in:\n "); 
    printf("S_E_T_U_P___T_R_E_N_D___L_I_N_E_A_R ( ) from Trend_Control.c\n");
    printf("The program will exit\n");
    exit(0);
  }
}

void S_E_T_U_P___T_R_E_N_D___M_U_L_T_I___S_T_E_P ( Trend_Control * T,
						   Parameter_Table * Table,
						   int Input_Parameter )
{
  double Value_0, Value_1;

  assert( No_of_TDC_FUNC_AUX_PARAM_MAX <= 10 );
  
  T->Tr_No_of_Jumps = (int)Table->TDC->TDC_Function_Params[Input_Parameter][9];
  
  T->Tr_Time_0     = Table->T->Time_0;
  T->Tr_value_0    = Table->TDC->TDC_Function_Params[Input_Parameter][0];
  
  T->Tr_Time_1     = Table->TDC->TDC_Function_Params[Input_Parameter][1];
  T->Tr_value_1    = Table->TDC->TDC_Function_Params[Input_Parameter][2];
  
  T->Tr_Time_2     = Table->TDC->TDC_Function_Params[Input_Parameter][3];
  T->Tr_value_2    = Table->TDC->TDC_Function_Params[Input_Parameter][4];
  
  T->Tr_Time_3     = Table->TDC->TDC_Function_Params[Input_Parameter][5];
  T->Tr_value_3    = Table->TDC->TDC_Function_Params[Input_Parameter][6];
  
  T->Tr_Time_4     = Table->TDC->TDC_Function_Params[Input_Parameter][7];
  T->Tr_value_4    = Table->TDC->TDC_Function_Params[Input_Parameter][8];

  T->Tr_Time_5     = Table->T->Time_1;

  T->Tr_Input_Parameter = Input_Parameter;

  T->TYPE_of_TREND  = 5;

  Table->Tr         = T;
}

void Upload_Argument_Input_Trend_Values_into_Table ( Parameter_Table * Table,
						     int Input_Parameter, 
						     int pattern )
{
  switch (pattern)
  {
  case 5: //(multi-step function)
      
    Table->TDC->TDC_Function_Params[Input_Parameter][0] = Tr_value_0;
  
    Table->TDC->TDC_Function_Params[Input_Parameter][1] = Tr_Time_1;
    Table->TDC->TDC_Function_Params[Input_Parameter][2] = Tr_value_1;
    
    Table->TDC->TDC_Function_Params[Input_Parameter][3] = Tr_Time_2;
    Table->TDC->TDC_Function_Params[Input_Parameter][4] = Tr_value_2;
    
    Table->TDC->TDC_Function_Params[Input_Parameter][5] = Tr_Time_3;
    Table->TDC->TDC_Function_Params[Input_Parameter][6] = Tr_value_3;
    
    Table->TDC->TDC_Function_Params[Input_Parameter][7] = Tr_Time_4;
    Table->TDC->TDC_Function_Params[Input_Parameter][8] = Tr_value_4;
    
    Table->TDC->TDC_Function_Params[Input_Parameter][9] = (double)Tr_No_of_Jumps;

  break;
  
  default :
    printf(" Trend_Control.c: Type of Trend is not defined\n");
    printf(" Trend pattern can be only 5\n");
    printf(" This function translates dummy parameters introduced\n");
    printf(" through command line into the member of Parameter_Table\n");
    printf(" that controls parametric time dependencies\n"); 
    printf(" but here Trend Type is %d\n", pattern); Print_Press_Key(1,0,".");
    exit(0);
  
  }
}

void  Upload_Auxiliary_Parameter_Values_into_Table ( Parameter_Table * Table,
						     int Input_Parameter, 
						     int pattern )
{
  /* The definition of default values for auxiliary/secondary parameters, this is, 
     the ones that control the time dependence of the 'MODEL_PARAMETERS_MAXIMUM' 
     true model paremteres through pattern functions is fully context dependent. 
     This means that the actual values to use will depend on which function drives this 
     dependence and which parameters control this function. An array of 
     'No_of_TDC_FUNC_AUX_PARAM_MAX' secondary parameters is associated to every 
     true model parameter. The actual values for No_of_TDC_FUNC_AUX_PARAM_MAX and 
     MODEL_PARAMETERS_MAXIMUM are defiend in MODEL.h 
  */
  
  int i;
  
  if (Input_Parameter == Tr_Input_Parameter ) {
    /* In this case, the values controling functional trend 
       are speficified though the comand line */
    Upload_Argument_Input_Trend_Values_into_Table ( Table,
						    Input_Parameter, 
						    pattern );
  }
  else {

    switch(Input_Parameter)
      {
      case  0: /* Mu_0 */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 	
           
	break;
      case  1: /* Gamma_0 */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.5;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 20.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 2.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 1520.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 1.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 2000.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 2.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 2100.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 1.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 5.0;

	printf(" Changing value in Gamma_0 trend\n");
	Print_Press_Key(1,0,".");
	
	break;
      case  2: /* k_0  */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
	break;
      case  3: /* Mu_1 */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
	break;
      case  4: /* Gamma_1 */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
	break;
      case  5: /* k_1 */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
	break;
      case  6: /* Mu_2 */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
	break;
      case  7: /* Gamma_2 */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
	break;
      case  8: /* k_2 */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
	break;
      case  9: /* Mu_3 */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
	break;
      case 10: /* Gamma_3 */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
	break;
      case 11: /* k_3 */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
      break;
      case 12: /* Mu_4 */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
	break;
      case 13: /* Gamma_4 */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
	break;
      case 14: /* k_4 */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
	break;
      case 15: /* No_of_GROUPS */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
	break;
      case 16: /* P->A_Rate */
	Table->TDC->TDC_Function_Params[Input_Parameter][0] = 0.0;
  
	Table->TDC->TDC_Function_Params[Input_Parameter][1] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][2] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][3] = 0.0;
	Table->TDC->TDC_Function_Params[Input_Parameter][4] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][5] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][6] = 0.0; 
	
	Table->TDC->TDC_Function_Params[Input_Parameter][7] = 0.0; 
	Table->TDC->TDC_Function_Params[Input_Parameter][8] = 0.0; 
    
	Table->TDC->TDC_Function_Params[Input_Parameter][9] = 0.0; 
	break;
	
      default:
	printf(".... INVALID PARAMETER KEY (key = %d)\n", Input_Parameter);

	printf(".... The permited correspondences are:\n");
	printf("\n");
	fprintf_Model_Parameters(stdout, Table);
	
	printf(" The maximum number of parameters is Number_PAR_MAX\n");
	printf(" The permited number of keys go from 0, to %d\n", MODEL_PARAMETERS_MAXIMUM-1);
	
	exit(0);
      }

    
  }
  
}
  
