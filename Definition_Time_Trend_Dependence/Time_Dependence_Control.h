void Time_Dependence_Apply_Optimized(Parameter_Table * Table, int j);

void Time_Dependence_Apply(Parameter_Table * , double );

int Integer_Position_of_a_Time(Parameter_Table * , double );

double Time_Dependence_Resolve(Parameter_Table * , int , int , double );

double Time_Dependence_Sigmoidal_Function ( Parameter_Table * , int , double ); 

double Time_Dependence_Sinusoidal_Function ( Parameter_Table * , int , double ); 

double  Triangular_Trend_Function (Parameter_Table * , double );

double  Linear_Trend_Function (Parameter_Table * , double );

double Multi_Step_Function (Parameter_Table * , double );

double Sinusoidal_Oscillation(double , double , double , double , double );

void Time_Dependence_Control_Alloc ( Time_Control * ,
				     Time_Dependence_Control * ,
				     Parameter_Table * ,
				     int ,
				     int ,
				     int );
void Time_Dependence_Control_Upload (Time_Control * ,
				     Time_Dependence_Control * ,
				     Parameter_Table * ,
				     int ,
				     int ,
				     int ,
				     int ,
				     int ,
				     int ,
				     int ,
				     int ,
				     int * ,
				     int * ,
				     char * ,
				     char * );

void Time_Dependence_Control_Upload_Optimized (Time_Control * ,
					       Time_Dependence_Control * ,
					       Parameter_Table * ,
					       int ,
					       int , 
					       int ,
					       int ,
					       int ,
					       int ,
					       int ,
					       int ,
					       int * ,
					       int * ,
					       double ** ,
					       double ** ,
					       double  * );

void Time_Dependence_Control_Free ( Time_Dependence_Control * , Parameter_Table * ); 
