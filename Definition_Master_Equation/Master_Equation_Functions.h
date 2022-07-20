int jacobian_ME (double t, const double y[], double *dfdy, double dfdt[], void *params);

int function_ME (double t, const double y[], double dydt[], void *params);

void JACOBIAN_Matrix_ME(gsl_matrix * m, const double *y, double t, int W, Parameter_Table * P);

void Model_Parameters_Master_Equation(Parameter_Table * Table,
				      int * No_of_CONFIGURATIONAL_STATES,
				      int * n_DIMENSION,
				      int * n_x, int * n_y, int * n_z); 

void i_to_nm_Map(Parameter_Table * Table,
		 int i, int * n_i, int * m_i);

void nm_to_i_Map(Parameter_Table * Table,
		 int * i, int n, int m);

double Sumando(int a_0, int m);

void Stationary_Probability_Distribution(Parameter_Table * Table);

void Stationary_Distribution_Allocation_Initialization( Parameter_Table * Table );

void Stationary_Distribution_Free ( Parameter_Table * Table ); 

double PS_nn_Function( int n, int m, double C0, double p );

double PS_nk_Function( int n, int k, double Cn, double q );

void Marginal_Stationary_Probabilities_Calculation ( Parameter_Table * Table );

void Probability_Distribution_Vector_into_Matrix_Form( Master_Equation * ME );

void Labels_for_Marginal_Probabilities (Parameter_Table * Table); 

void Marginal_Probability_Calculation ( Parameter_Table * Table );

void Marginal_Probability_Averages_Calculation ( Parameter_Table * Table );

void Print_Probability_Distribution ( Parameter_Table * Table ); 

void Print_Marginal_Averages( double Time_Current, Parameter_Table * Table); 

int master_equation_driver( Parameter_Table * Table, int j, double * Time_Current );

void Normalization_Master_Equation(Parameter_Table * Table); 

int master_equation_time_dynamics( Parameter_Table * Table );

void Initial_Condition_Master_Equation( Parameter_Table * Table, double * y_INI );

void Common_Initial_Condition_Command_Line_Arguments_into_Table(Parameter_Table * Table); 

void C_P_G___M_A_R_G_I_N_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Parameter_Table * Table,
							 int j,
							 int n, 
							 double Time_Current,
							 int SAME );

void C_P_G___E_M_P_I_R_I_C_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Parameter_Table * Table,
							   int j,
							   int n, 
							   double Time_Current,
							   int SAME );

void C_P_G___S_T_A_T_I_O_N_A_R_Y___D_I_S_T_R_I_B_U_T_I_O_N ( Parameter_Table * Table, 
							     int n, 
							     int SAME );

void Saving_Marginal_Distribution(Parameter_Table * Table, int j, int n, double Time_Current); 

void Saving_Marginal_Distribution_Triplets(Parameter_Table * Table, int j, double Time_Current); 

void Saving_Empirical_Distribution_vs_ME_Numerical_Integration ( Parameter_Table * Table,
								 int j,
								 double Time_Current);
