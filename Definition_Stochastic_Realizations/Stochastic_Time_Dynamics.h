int S_T_O_C_H_A_S_T_I_C___T_I_M_E___D_Y_N_A_M_I_C_S( int i, 
						      Parameter_Table * Table, 
						      int * Bad_Times );

int Stochastic_Time_Dynamics_Numerical( int i,
					Parameter_Table * Table,
					int * Bad_Times );

void Initial_Conditions_Stochastic_Dynamics( Parameter_Table * Table, double * y_INI );

void Get_Initial_y_INI_Random_Vector_into_Integer_Values(Parameter_Table * Table,
							 double * y_INI); 

void Initial_Conditions_Global_Populations (Parameter_Table * Table, double * Vector_Model_Variables );

double Bacterial_Type_Population_per_Cell (Parameter_Table * Table, int i_Strain, int j_Patch);

double Plasmid_Type_Population_per_Cell (Parameter_Table * Table, int i_Plasmid, int j_Patch);

void Patch_System_Initialization (Community ** PATCH, Parameter_Table * Table, double * y_INI);

int Advance_Current_Time( Parameter_Table * Table, 
			   Stochastic_Rate * Rate, double * Time_Current, int * New );

int Choose_Village(double max_Probability, Community ** Pop, Parameter_Model * Par);

int Choose_Village_Binary_Tree(double max_Probability, Community ** Pop, Parameter_Model * Par);

void Choose_Village_and_Event_Binary_Tree( double max_Probability, Community ** Pop,
                                           Parameter_Model * Par, 
                                           int * patch, int * event );

void Choose_Village_and_Event_Priority_Queu(double max_Probability, Community ** Pop,
                                       Parameter_Model * Par, 
                                       int * patch, int * event );

void Execute_One_Step(Community ** SP,
		      Parameter_Table * Table,
		      double max_Probability, 
		      int * Event, int * Patch);

void Temporal_Dynamics(Community ** My_Community, Parameter_Table * Table, Stochastic_Rate * Rate);

void Temporal_Dynamics_Update( Community ** My_Community,
			       Parameter_Table * Table, Stochastic_Rate * Rate,
			       int Type_of_Event, int * Patch);

void Local_Population_Decrease ( int n, Community * Patch );

void Local_Population_Increase ( int n, Community * Patch );

void Some_Other_Patch_Population_Decrease(int , int , Parameter_Table * );

int Some_Other_Patch_Population_Increase(int , int , Parameter_Table * );

int Some_Other_Species_Establishes(int , int S, Parameter_Table * Table);

int Other_Patch_Population_Increase(int Event, int x, int n, int Sp,
					                Parameter_Table * Table);

int Boundness_Control( int Event, Parameter_Table * Table,
			                 int x_out, int x_in, int Sp, int jS, double Y, int J);

void Boundness_Control_Local_Patch(int Event, Parameter_Table * Table, 
                                   int x, int jS, double Y, int J);							 

void Positivity_Control( int Event, Parameter_Table * Table,
			 int x, int jS, double Y, int J);

void Event_Adjacence_List_Initialization(Community ** PATCH,
					 Parameter_Model * P); 

void Event_Delta_Matrix_Initialization(Community ** PATCH,
				       Parameter_Model * P); 

