void Write_Parameter_Table ( Parameter_Table * P, int no_Par ); 

void write_Parameter_Table___RANGES___VALUES___LATEX ( char * File_Model_Parameters,
						       char * Type_of_Boundary_Parameter_Space,
						       Parameter_Table * P,
						       double * Par_MAX, double * Par_min,
						       int no_Par ); 

void write_Param_Table_PARAMETER_RANGES___LATEX ( char * File_Model_Parameters,
						  char * Type_of_Boundary_Parameter_Space,
						  Parameter_Table * P, 
						  double * Par_MAX, double * Par_min,
						  int no_Par );

void write_Param_Table_AVE_STD___LATEX( char * File_Model_Parameters,
					Parameter_Table * P, 
					double * Par_AVE, double * Par_STD, int no_Par ); 

void write_Param_Table_AVE_CI___LATEX( char * File_Model_Parameters,
				       Parameter_Table * P, 
				       double * Par_AVE,
				       double ** Confidence_Interval,
				       int no_Par );
