void Parameter_Space_Boundary_for_Trends(int j, double *value_0, double *value_1);

void Initial_Condition_Boundary_Values(int j, double *value_0, double *value_1, Parameter_Table * P);

#if defined PARAMETRIC_CONFIGURATIONS
void Parameter_Space_Alloc( Parameter_Space * , int , int , int *  );
#else
void Parameter_Space_Alloc( Parameter_Space * , int , int *  );
#endif

void Parameter_Space_Initialization(Parameter_Space * , int , double, int,
				    int * , int * , double ** , double * );

void Parameter_Space_IC_Initialization(Parameter_Space * , int ,
				       int * , int * , double ** , double * );

void Parameter_Space_Error_Initialization(Parameter_Space * , int ,
				       int * , int * , double ** , double * );


void Parameter_Space_Free(Parameter_Space *, int );


void Parameter_Index_Checking_Ordering (int * , int * , int ,
					int * , int * , int * , int * ,
					int * , int * , int * );

void Parameter_Space_Write_Min_Max_IC_Values (Parameter_Space * , Parameter_Table *  );

void Parameter_Space_Write_Min_Max_Values (Parameter_Space * , Parameter_Table *  );

void Parameter_Space_Write_Min_Max_Error_Values (Parameter_Space * , Parameter_Table * ); 


void Random_Initial_Guess_within_Boundaries_Model (Parameter_Model * Initial_Guess,
						   Parameter_Space * S );

void Random_Initial_Guess_within_Boundaries_Table (Parameter_Table * Initial_Guess,
						   Parameter_Space * S ); 

/* The next two functions are not generic. They cannot be transportable as library functions: they are genuine of the AIDS project */
int Random_Initial_Guess_within_Boundaries_and_Data (Parameter_Table * Initial_Guess, 
						     Parameter_Space * S, double ** Data,
						     int Alpha_Index );

double Alpha_from_Sigma_Parameters ( double x_e, Parameter_Table * Table ); 


