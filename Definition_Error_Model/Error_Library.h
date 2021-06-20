typedef struct Parameter_Table_Root_Solverinfo
{

  gsl_function * f;
  double      r;

}Parameter_Table_Root_Solver;

void my_error_handler (const char * reason, const char * file, int line, int gsl_errno); 

double Measurement_Error_Function (gsl_rng * r, Parameter_Table * P,  double y);

double my_Density_Function (double x, void * params);

double Cummulative_Density_Function( double x_lo, double x_hi, gsl_function * Density_Function );

double Cummulative_Distribution_Function( double x, gsl_function * Density_Function );

double Function_Root_Solver( double x, void * p );

double da_gsl_ran_continuous_Function( const gsl_rng * r, gsl_function * Density_Function, int i, int N );

void Parameter_Table_into_Vector_Entries_Error_Model ( Parameter_Table * P, gsl_vector * X,
						       int * Index,
						       int No_of_PARAMETERS,
						       int No_of_IC,
						       int No_of_Error_Parameters );

double Error_Model_into_Vector_Entry_Table(int j, Parameter_Table * Table);


double Error_Model_into_Vector_Entry_Model(int j, Parameter_Model * Table); 


void Vector_Entries_into_Parameter_Table_Error_Model ( const gsl_vector * X,
						       Parameter_Table * P,
						       int * Index,
						       int No_of_PARAMETERS,
						       int No_of_IC,
						       int No_of_Error_Parameters);

void Vector_Entry_into_Error_Model_Table(double value, int j, Parameter_Table * Table);

void Vector_Entry_into_Error_Model_Model(double value, int j, Parameter_Model * Table); 

void Random_Error_Model_within_Boundaries_Model (Parameter_Model * ,
						 Parameter_Space *  ); 

void Random_Error_Model_within_Boundaries_Table (Parameter_Table * ,
						 Parameter_Space *  ); 
