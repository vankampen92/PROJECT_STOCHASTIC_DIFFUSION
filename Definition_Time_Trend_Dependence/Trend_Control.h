void T_R_E_N_D___C_O_N_T_R_O_L___F_R_E_E( Trend_Control * Tr );

void T_R_E_N_D___C_O_N_T_R_O_L___A_L_L_O_C( Trend_Control * Tr, Parameter_Table * P );

void  T_R_E_N_D___C_O_N_T_R_O_L___U_P_L_O_A_D( Trend_Control * T, Parameter_Table * Table );


void S_E_T_U_P___T_R_E_N_D___T_R_I_A_N_G_L_E ( Trend_Control * T, 
					       Time_Control * Time, 
					       int Input_Parameter );

void S_E_T_U_P___T_R_E_N_D___L_I_N_E_A_R ( Trend_Control * T, 
					   Time_Control * Time, 
					   int Input_Parameter, 
					   int Signature );

void S_E_T_U_P___T_R_E_N_D___M_U_L_T_I___S_T_E_P ( Trend_Control * T,
						   Parameter_Table * Table,
						   int Input_Parameter ); 

void Upload_Argument_Input_Trend_Values_into_Table ( Parameter_Table * Table,
						     int Input_Parameter, 
						     int pattern ); 

void  Upload_Auxiliary_Parameter_Values_into_Table ( Parameter_Table * Table,
						     int Input_Parameter, 
						     int pattern ); 
