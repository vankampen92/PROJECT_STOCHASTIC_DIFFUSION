void  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N (Parameter_Table * Table, Parameter_Model * P);

void P_A_R_A_M_E_T_E_R___M_O_D_E_L___F_R_E_E( Parameter_Model * P);

void Parameter_Model_Copy (Parameter_Model * P_Destination, Parameter_Model * P_Source);

void Vector_Entries_into_Parameter_Model ( const gsl_vector * X, Parameter_Model * P,
					   int * Parameter_Index, int No_of_PARAMETERS );

void Vector_Entry_into_Parameter_Model ( double value,
					 int key, Parameter_Model * P );

void Parameter_Model_into_Vector_Entries ( Parameter_Model * P, gsl_vector * X,
					   int * Parameter_Index, int No_of_PARAMETERS );

double Parameter_Model_into_Vector_Entry ( int key, Parameter_Model * P );

void Parameter_Model_Copy_into_Parameter_Table (Parameter_Table * P_Destination, Parameter_Model * P_Source);

void Parameter_Table_Copy_into_Parameter_Model (Parameter_Model * P_Destination, Parameter_Table * P_Source);
