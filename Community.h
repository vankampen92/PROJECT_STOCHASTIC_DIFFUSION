void Community_Allocation ( Community ** , Parameter_Model *  );
void Community_Free (Community ** , Parameter_Model * );
void Community_Initialization (Community ** , Parameter_Model * );
void Immigration_Preassure_on_Focal_Patch_Initialization( Community ** , Parameter_Model * P ); 
void Network_Structure_Inititialization (Community ** , int, int );
void Writing_Adjacency_List(Community ** );
void Print_Meta_Community_Patch_System (Parameter_Table *);
void Writing_Adjacency_List_VonNeumann(Community ** PATCH);
void Set_Von_Neumann_1st_Neighbors(Community ** PATCH, int no, int N_X, int N_Y, int i); 
void Community_Scatter_Plot_Representation( Parameter_Table * Table,
					    int i_Replicate, int j_Time);
void Community_Binary_Tree_Initialization (Parameter_Table * Table);
