int generic_Function_Parameter_2Dim_Scan( Parameter_Table * P,
                                          int No_of_POINTS_1, int Input_Parameter_1,
                                          int No_of_POINTS_2, int Input_Parameter_2,
                                          double (* GENERIC_FUNCTION) (Parameter_Table *),
                                          double * W_GRID,
                                          char * generic_Function_Parameter_Scan_Output_File );

int generic_Function_Parameter_2Dim_Scan_Improved( Parameter_Table * P,
                                          int No_of_POINTS_1, int Input_Parameter_1,
                                          int No_of_POINTS_2, int Input_Parameter_2,
                                          double (* GENERIC_FUNCTION) (Parameter_Table *),
                                          double * W_GRID,
                                          char * generic_Function_Parameter_Scan_Output_File,
					  int X_LINEAR, int Y_LINEAR);

double Function_to_Minimize( Parameter_Table * );

void Confidence_Intervals_from_Likelihood_Profile( double * Profile, double * X, int N,
						   double LIKELIHOOD_JUMP, 
						   double * X_MLE, double * CI ); 

double GSL_Function_to_Minimize( const gsl_vector * x, void * Par );

double GSL_Function_to_Minimize_AIDS( const gsl_vector * x, void * Par );

double GSL_Function_to_Minimize_Error_Model( const gsl_vector * x, void * Par );

double GSL_Function_to_Minimize_Binomial_Free_Consumers( const gsl_vector * x, void * Par );

double GSL_Function_to_Minimize_Multinomial_Free_Consumers( const gsl_vector * x, void * Par );

void Theta_Nu_Sums (double * Theta_Vector, double * Nu_Vector, 
                    double * Theta_Sum, double * Sum_ThetaNu_Ratio, 
                    int N);

double Function___poft (double Theta, double Nu, 
                        double Sum_ThetaNu_Ratio, double Theta_Sum, 
                        double t);

double GSL_Function_to_Minimize_Beddington_DeAngelis( const gsl_vector * x, void * Par ); 

double GSL_Function_to_Minimize_Beddington_DeAngelis_Marginal_0(const gsl_vector * x, void * Par);

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

/* Eigenvalue Calculation */
void GSL_Eigenvalue_Calculation ( double * y_Sol,  int N,
                                  Parameter_Table * P,
                                  double * l_re, double * l_im );

void Dominant_Eigenvalue_Calculation(double * Y1, double * Y2, int N,
				     int * Index_Value_D, int * Index_Value_S);

void NR_Eigenvalue_Calculation ( double * y_Sol,  int N, 
				 Parameter_Table * P, 
				 double * l_re, double * l_im ); 

void NR_Eigenvalue_Calculation_float ( double * y_Sol,  int K, int W, 
				       Parameter_Table * P, 
				       float * l_re, float * l_im ); 

double Well_Defined_Matrix_Elements(float **mm, int N, float xmin, float xmax);

int showing_eigenValues(float *l_re, float *l_im, int n); 

int fprintf_to_File_Matrix_gsl(FILE * Out, gsl_matrix * A, int MM, int NN); 

int gsl_matrix_to_NR_matrix(gsl_matrix * A, float **a,  int MM, int NN);

/* Backup Version */
void E_I_G_E_N___V_A_L_U_E___C_A_L_C_U_L_A_T_I_O_N ( double * y_Sol,  int K, int W,
                                                     Parameter_Table * P,
                                                     float * l_re, float * l_im );
#include "treenode.h"

