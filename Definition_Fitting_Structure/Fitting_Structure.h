void Parameter_Fitting_Alloc( Parameter_Fitting * F, int , Parameter_Table * );

void Parameter_Fitting_Initialization( Parameter_Fitting * F, int ,
				       Observed_Data * , Parameter_Table * ); 

void Parameter_Fitting_Free(Parameter_Fitting * F); 

void Parametric_Configurations_into_Fitting_Structure_from_File (Parameter_Fitting * ,
								 char *  ); 

void Parametric_Configurations_from_Fitting_Structure_into_File (Parameter_Fitting * ,
								 char *  ,
								 int ); 
