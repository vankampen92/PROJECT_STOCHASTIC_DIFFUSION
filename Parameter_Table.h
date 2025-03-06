void P_A_R_A_M_E_T_E_R___T_A_B_L_E___A_L_L_O_C( Parameter_Table * );

void P_A_R_A_M_E_T_E_R___T_A_B_L_E___U_P_L_O_A_D( Parameter_Table *, int * );

void P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( Parameter_Table * );

void Parameter_Table_Index_Update(int * , int , Parameter_Table * );

void Parameter_Values_into_Parameter_Table(Parameter_Table * ); 

void Resetting_Lambda_Delta_Vectors (Parameter_Table * ); 

void Resetting_Alpha_Nu_Vectors (Parameter_Table * ); 

void Resetting_Alpha_Nu_Vectors_Constant (Parameter_Table * Table);

void Resetting_Multiresource_Levels (Parameter_Table * Table);

void Writing_Alpha_Nu_Theta_Vectors(Parameter_Table * Table);

void Ressetting_Species_Characteristic_Parameters (Parameter_Table * Table);

void Setting_Interaction_Matrices (Parameter_Table * Table);

int Determining_actual_No_of_RESOURCES(Parameter_Table * Parameter_Table);

void Setting_Adjacency_Lists_from_Interaction_Matrices (Parameter_Table *Table);

void Calculate_Strain_and_Profile(Parameter_Table * Table, int n, int * i_Strain, int * k_Profile);

void Printing_Strains_Profiles_and_Lists(Parameter_Table * Table);

void Setting_Plasmid_Characteristic_Parameters (Parameter_Table * Table);

void Setting_Strain_Characteristic_Parameters (Parameter_Table * Table);

void Setting_up_Constant_Metapopulation_Connectivity_Matrix (Parameter_Table * Table);
