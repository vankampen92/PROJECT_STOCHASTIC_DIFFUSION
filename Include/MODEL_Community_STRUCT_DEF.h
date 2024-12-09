typedef struct Communityinfo
{
  /* Total population of each type of class */
  int No_of_RESOURCES;          
    
  int No_of_CELLS;   /* Number of patches of the system the single patch is embedded in         */
  int No_of_CELLS_X; /* Number of patches of the system the single patch is embedded in (X Dim) */
  int No_of_CELLS_Y; /* Number of patches of the system the single patch is embedded in (Y Dim  */
  double X_DIMENSION;
  double Y_DIMENSION; 
  
  int LOCAL_STATE_VARIABLES; /* Number of variables needed to describe the state
			                          of the patch at any given time. This only coicides with the 
			                          total MODEL_STATE_VARIABLES when the we have only one patch     */

  int m_0;   /* Total Empty Space in a local community                                          */
  int * n;   /* n[0], ..., n[LOCAL_STATE_VARIABLES-1] defines the state of the patch completely */
  double * rate;     /* Transition probability for just an individual of each type         */
  double * rToI;     /* Transition probability for any among all individuals of each type  */ 
  double ratePatch;  /* Transition probability of this patch                               */
  
  struct point center;  /* Coordinates of the position of the center of patch    */

  int No_NEI;

  double ** Out_Migration_Vector;
  double ** In_Migration_Vector;
	 
  double * Total_Per_Capita_Out_Migration_Rate;
  double * Total_Imm_Rate_Preassure;

  double ** Imm_Rates_Preassure;
      
  int * Patch_Connections; 

  double *** Metapop_Connectivity_Matrix; 

  int ** Event_Adjacence_List;

  double ** Event_Delta_Matrix;

  double *** Event_Delta_Tensor;  /* One matrix per species */

  struct Communityinfo ** NEI; /* An array of pointers to the neighbors the focal patch is connected to */

  Strain ** Local_Strain_Population;     /* Differents strains subpopulations thriving in a local community 
                                           (DIFFUSION_ECO_PLASMIDS) 
                                         */
  Plasmid ** Local_Plasmid_Population;   /* (DIFFUSION_ECO_PLASMIDS)
                                         */
  int * Bacterial_Type_Population; /* Local total number of individual cells per bacterial type (regardless profile) */

  int * Plasmid_Type_Population;   /* Local number of plasmid types (across all different bacterial types) */  
}Community;

