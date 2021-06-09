int generic_Function_Parameter_2Dim_Scan( Parameter_Table * P,
                                          int No_of_POINTS_1, int Input_Parameter_1,
                                          int No_of_POINTS_2, int Input_Parameter_2,
                                          double (* GENERIC_FUNCTION) (Parameter_Table *),
                                          double * W_GRID,
                                          char * generic_Function_Parameter_Scan_Output_File );

double Function_to_Minimize( Parameter_Table * );


double GSL_Function_to_Minimize( const gsl_vector * x, void * Par );

double GSL_Function_to_Minimize_AIDS( const gsl_vector * x, void * Par );

double GSL_Function_to_Minimize_Error_Model( const gsl_vector * x, void * Par );

double GSL_neglog_Error_Probability_Model( double * Data, double * Theory,
					   int N , int No_of_VARIABLES,
					   Parameter_Fitting  * F, 
	                                   double (* Error_Model )(double, double, Parameter_Fitting *) ); 

double GSL_neglog_Error_Probability_Model_Gaussian( double Data, double Theory,
						    Parameter_Fitting  * F );

double GSL_neglog_Error_Probability_Model_Gamma( double Data, double Theory,
						 Parameter_Fitting  * F ); 

double GSL_Minimization_Driver( Parameter_Fitting * F );

double GSL_Minimization_Simplex (Parameter_Fitting * F,
				 gsl_vector * Initial_Guess,
				 gsl_vector * Solution,
				 double ( * Function )( const gsl_vector * , void * ) );

void GSL_CPGPLOT_Minimization_Simplex (Parameter_Fitting * F,
				       gsl_vector * Solution, size_t iter,
				       double ( * Function )( const gsl_vector * , void * ) );

int Checking_for_Parameter_Boundaries( Parameter_Fitting * F, const gsl_vector * x );

double Inspecting_Likelihood_of_Final_Solution( const gsl_vector * x, void * Par ); 

double Inspecting_Solution_Driver( Parameter_Fitting * F ); 
