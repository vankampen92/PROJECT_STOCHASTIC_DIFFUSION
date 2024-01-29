void Fixed_Points_All( Parameter_Table * Table,
		       double * Vector_Stationarity_Lower,
		       double * Vector_Stationarity_Inter,
		       double * Vector_Stationarity_Upper,
		       double Epsilon);

int  Coexistence_Condition ( Parameter_Table * Table );

double Coexistence_Condition_Double ( Parameter_Table * Table ); 

double Function_to_Type_of_Stability( Parameter_Table * Table ); 

void Stationary_Solution_Feasibility_Control (Parameter_Table * Table ); 

void Stationary_Solution_dvdt_Check (Parameter_Table * Table );

void assert_positive_model_parameters (Parameter_Table * Table ); 

void Fixed_Points_Linear_System_HII_nD(Parameter_Table * Table, gsl_vector * y);

double R0_Function_W( Parameter_Table * Table );

double R0_Function_F( Parameter_Table * Table );