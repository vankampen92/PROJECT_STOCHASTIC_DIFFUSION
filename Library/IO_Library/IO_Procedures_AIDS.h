void Reading_Demographic_Parameters_from_File(char * , double ** , int * , int ); 
   
void Reading_Model_Parameters_from_File(char * , double ** , int * , int ); 

void Writing_Demographic_Parameters_Matrix(Parameter_Table * , const char * , double ** , int ,
					   int * , int );

void Writing_Model_Parameters_Matrix(Parameter_Table * , const char * ,
				     double ** ,
				     int , int ); 

void Uploading_Demographic_Parameters_into_Parameter_Table(Parameter_Table * , double ** , int ,
							   int * , int  ); 

void Uploading_Model_Parameters_into_Parameter_Table(Parameter_Table * , double ** , int k); 

void Reading_Observed_Data(char * , double ** , int , int , int * , int * , int,
			   double * , double * , double *);
