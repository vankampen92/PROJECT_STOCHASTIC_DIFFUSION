void Reading_Model_Parameters_from_File(char * , double ** , int * , int ); 

void Writing_Model_Parameters_into_File(char * , char ** , double ** , int ,
					int ); 

void Writing_Model_Parameters_Matrix(Parameter_Table * , const char * ,
				     double ** ,
				     int , int ); 

void Creating_Title_Row (Parameter_Table * Table, char ** Title_Parameters); 

void Uploading_Model_Parameters_into_Parameter_Table(Parameter_Table * , double ** , int k); 

void Uploading_Full_Model_Parameters_into_Parameter_Table(Parameter_Table * Table, double ** Data,
							  int k); 

void Reading_Observed_Data(char * , double ** , int , int , int * , int * , int,
			   double * , double * , double *);

