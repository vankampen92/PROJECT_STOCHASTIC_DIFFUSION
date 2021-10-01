int jacobian_ME (double t, const double y[], double *dfdy, double dfdt[], void *params);

int function_ME (double t, const double y[], double dydt[], void *params);

void JACOBIAN_Matrix_ME(gsl_matrix * m, const double *y, double t, int W, Parameter_Table * P);

void Model_Parameters_Master_Equation(Parameter_Table * Table,
				      int * No_of_CONFIGURATIONAL_STATES,
				      int * n_DIMENSION,
				      int * n_x, int * n_y, int * n_z); 

double Average_Number_of_Feeding_Events( Parameter_Table * Table, double T ); 

int master_equation_driver( Parameter_Table * Table, int j, double * Time_Current );

void Normalization_Master_Equation(Parameter_Table * Table); 

int master_equation_time_dynamics( Parameter_Table * Table );

void Initial_Condition_Master_Equation( Parameter_Table * Table, double * y_INI );


