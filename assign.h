void Vector_Entries_into_Parameter_Table ( const gsl_vector * X, Parameter_Table * P,
					   int * Parameter_Index, int No_of_PARAMETERS );

void Vector_Entries_into_Parameter_Table_Initial_Condition ( const gsl_vector * X,
							     Parameter_Table * P,
							     int * Index,
							     int No_of_PARAMETERS,
							     int No_of_IC); 

void AssignVectorEntry_to_Structure(Parameter_Table * P, int j, double value);

void Parameter_Table_into_Vector_Entries ( Parameter_Table * P, gsl_vector * X,
					   int * Parameter_Index, int No_of_PARAMETERS );

void Parameter_Table_into_Vector_Entries_Initial_Condition ( Parameter_Table * P, gsl_vector * X,
							     int * Index,
							     int No_of_PARAMETERS,
							     int No_of_IC); 

double AssignStructValue_to_VectorEntry(int j, Parameter_Table * P);

void AssignCodes_to_Model_Parameters(int j, char * Label, Parameter_Table *P);

void AssignLabel_to_Model_Parameters(int j, char * Label, Parameter_Table *P);
void Label_to_Model_Parameters(int j, char * Label, Parameter_Table *P);

void AssignSymbol_to_Model_Parameters(int j, char * Label, Parameter_Table *P);
void Symbol_to_Model_Parameters(int j, char * Label, Parameter_Table * P); 

void AssignCPGPLOT_Symbol_to_Model_Parameters(int j, char * Label, Parameter_Table *P);
void CPGPLOT_Symbol_to_Model_Parameters(int j, char * Label, Parameter_Table * P); 

void AssignLabel_to_Model_Parameters__LATEX(int j, char * Label, Parameter_Table *P);
void Label_to_Model_Parameters__LATEX(int j, char * Label, Parameter_Table *P);

void AssignLabel_to_Model_Parameters__LATEX__SYMBOL(int j, char * Label, Parameter_Table *P);
void Label_to_Model_Parameters__LATEX__SYMBOL(int j, char * Label, Parameter_Table *P); 

void AssignLabel_to_Model_Variables(int j, char * Label, Parameter_Table *P);
void Calculate_Strain_and_Profile(Parameter_Table * Table, int n, int * i_Strain, int * k_Profile); /* DIFFUSION_ECO_PLASMIDS */

void AssignLabel_to_Output_Variables(int j, char * Label, Parameter_Table * P);

void AssignLongLabel_to_Output_Variables(int j, char * Label, Parameter_Table * P); 

void AssignCPGPLOT_Symbol_to_Output_Variables(int j, char * Label, Parameter_Table * P);

void AssignSymbol_to_Model_Variables(int j, char * Label, Parameter_Table *P);

void AssignSymbol_to_Output_Variables(int j, char * Label, Parameter_Table * P);

void Vector_Entry_into_Model_Variable_Initial_Condition_Table(double Value, 
								int j, Parameter_Table * P);

void Vector_Entry_into_Model_Variable_Initial_Condition_Model(double Value, 
								int j, Parameter_Model * P);

double Model_Variable_Initial_Condition_into_Vector_Entry_Table(int j, Parameter_Table * P);

double Model_Variable_Initial_Condition_into_Vector_Entry_Model(int j, Parameter_Model * P);

