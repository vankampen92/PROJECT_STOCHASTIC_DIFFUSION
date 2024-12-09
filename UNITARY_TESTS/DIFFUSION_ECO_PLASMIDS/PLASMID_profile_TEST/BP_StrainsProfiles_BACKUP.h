typedef struct Donor_Recipient_Pairinfo
{
  int ne;     /* Strain ID identifier */ 
  double ra;  /* Rate                 */
  int N;   /* Total Number of Effective Donors */
}Donor_Recipient_Pair;

typedef struct Time_Controlinfo
{
  int I_Time;
  double Time_0;
  double Time_1;
  double Delta_T; 

  double * Time_Vector;
}Time_Control;

typedef struct Parameter_Tableinfo
{
  Time_Control * T;   

  int TYPE_of_INITIAL_CONDITION; 

  int LOCAL_STATE_VARIABLES; 
  int MODEL_STATE_VARIABLES;

  double * Vector_Model_Variables;

  /* Interaction Matrices (DIFFUSION_ECO_PLASMIDS)                     */        
  double ** ABB;      /* Strain-Strain Competition Matrix              */
  double ** HBB;      /* Strain-Strain Conjugation Matrix              */
  double ** IBP;      /* Interaction Infection Bacteria-Plasmid        */
  double ** CPP;      /* Plasmid-Plasmid Compatibility Matrix          */
  /* Adjancency lists corresponding to these matrices                  */
  int ** Competition_List_Indeces; 
  int ** Conjugation_List_Indeces;     
  int ** Recipient_List_Indeces;
  int ** Putative_Recipient_List_Indeces; 
  int ** Donor_List_Indeces;
  int ** Plasmid_Compatibility_Indeces; 

  /* No of actual viable profiles per strain. It can be different per each strain */
  int * n;                  /* (DIFFUSION_ECO_PLASMIDS) */
  int * n_0;                /* n_0[i]: Strain ID of the plasmid-free profile 
                                       associated to the i-th bacterial type */
  int * n_R;                /* n_R[i]: No of Recipients of Strain ID 'i' */

  int No_of_PLASMIDS; /* Used in DIFFUSION_ECO_PLASMIDS */

  int No_of_PROFILES; /* Used in DIFFUSION_ECO_PLASMIDS */
  
  int No_of_STRAINS;  /* Used in DIFFUSION_ECO_PLASMIDS */

  int No_of_RESOURCES;

  int No_of_CONJUGATION_EVENTS; 
  
  int K_R;             /* System Size (or Carrying Capacity)                 */
  double p_2;          /* Sparsity Parameter (a Matrix Connectance)                */ 
  double Beta_R;       /* CELL_DIVISION_RATE         */ 
  double p_1;          /* Segregation Error */  
  double Delta_R_0;    /* BASAL_DEATH_RATE           */ 
  double Delta_R_1;    /* STRESS_INDUCED_DEATH       */ 
  double Delta_C_0;    /* Competition-Induced Mortality (required to define the CBB matrix) */
  double Lambda_R_1;   /* Common Conjugation Donor-Recipient Encounter Rate  */
  double Mu;           /* DIFFUSION_RATE             */
  
  double Chi_C_0;      /* Plasmid Transmission Probability */
  double Alpha_C_0;    /* Plasmid reproduction cost  */
  double Nu_C_0;       /* Plasmid resistance         */      
  
  double * Alpha_C;   /* Plasmid Reproduction Costs (DIFFUSION_ECO_PLASMIDS)                  */
                      
  double * Nu_C;      /* Plasmid Resistance (DIFFUSION_ECO_PLASMIDS)                          */
                     
  double * Eta_RP;    /* Strain-Specific Conjugation/Encounter Rates (DIFFUSION_ECO_PLASMIDS) */

  double * Beta_AP;   /* Cell Division Rate of Single Bacterial Cell (DIFFUSION_ECO_PLASMIDS) */
  
  double * Delta_AP;  /* Per Capita Death Rates (DIFFUSION_ECO_PLASMIDS) */
  
  double * Mu_RP;     /* Bacterial Cell Diffusion/Movement Rates (DIFFUSION_ECO_PLASMIDS)     */   
                      
  double * Segregation_Error; /* Strain specific segregation errors (DIFFUSION_ECO_PLASMIDS) */

  int *** Strain_Profiles;              /* (DIFFUSION_ECO_PLASMIDS) */
  int ** StrainType_and_Profile;        /* [Strain Type][Profile No Specifier] */
  double ** Competition_Induced_Death;  /* Used in DIFFUSION_ECO_PLASMIDS */

  double * Conjugation_Gain_Rate;    /* Used in DIFFUSION_ECO_PLASMIDS */
  double * Conjugation_Loss_Rata;    /* Used in DIFFUSION_ECO_PLASMIDS */    

  Donor_Recipient_Pair **** DoRe;    /* Used in DIFFUSION_ECO_PLASMIDS */

  #if defined CPGPLOT_REPRESENTATION
  #include <include.CPG.global.h>
    Parameter_CPGPLOT * CPG;
    Parameter_CPGPLOT * CPG_STO; 
  #endif
}Parameter_Table;

void Setting_Interaction_Matrices (Parameter_Table * Table);
int Determining_actual_No_of_RESOURCES(Parameter_Table * Parameter_Table);
int SumUP_Profile ( int * Profile, int N );
void Setting_Adjacency_Lists_from_Interaction_Matrices (Parameter_Table *Table);
void Setting_Putatitive_Recipient_Lists_of_Potential_Trasconjugants (Parameter_Table * Table);
void Calculate_Strain_and_Profile(Parameter_Table * Table, int n, int * i_Strain, int * k_Profile);
void Printing_Strains_Profiles_and_Lists(Parameter_Table * Table);
void Printing_Putative_Recipient_Lists(Parameter_Table * Table);
void Printing_Strains_Profiles(Parameter_Table * Table);
void Print_Strain_Profile(int * Profile, int No_of_PLASMIDS);
void Printf_Infection_Profile(int * Profile, int i_List, Parameter_Table * Table);
void Create_Binary_Combination( int ** Binary_Combination, int N, int LENGTH );
void int_buffer_rec(int ** Number_List, int N, int * number, int n, int length);
void GSL_Init_Random_Seed(const gsl_rng * );
void GSL_Init_Random_Seed_from_File(const gsl_rng * );
void show_DoubleMatrix(double **M, int Nx, int Ny);
int Calculate_Strain_ID_from_Profile_and_Strain_Type(Parameter_Table * Table, int Sp_Strain, int * Profile );
int Compatibility_Plasmid_Sets_Assert( Parameter_Table * Table, 
                                       int * Plasmid_Set_0, int N_0, int * Plasmid_Set_1, int N_1 );
int Compatibility_Plasmid_Sets_Count( Parameter_Table * Table, 
                                      int * Plasmid_Set_0, int N_0, int * Plasmid_Set_1, int N_1 );
int Compatibility_Profiles_Assert (Parameter_Table * Table, int * Profile_0, int * Profile_1, int N);
int Infection_Condition_Assert(int i_Sp, int * Profile, Parameter_Table * Table);
int Profile_Selfconsistency_Assert(int i_Strain_Sp, int * Profile, Parameter_Table * Table);
int Potential_Donor_Recipient_Pair_Assert(Parameter_Table * Table, 
                                          int Donor_ID, int Donor_Sp, int k_D,
                                          int Recip_ID, int Recip_Sp, int k_R);
int Effective_Donor_Recipient_Reaction_Assert(Parameter_Table * Table, 
                                              int Donor_ID, int Donor_Sp, int * Donor_Set, int n_D, 
                                              int Recip_ID, int Recip_Sp, int k_R);
void Setting_Reactive_Recipient_Donor_Pairs_and_Rates(Parameter_Table * Table); 
void Transconjugation_Gain_and_Loss_Total_Rates( double * Y, double * Gain, double * Loss, 
                                                 Parameter_Table * Table );               /* In use when defining the ODEs */
int Recipient_Donor_Transconjugant_Rate( int Trans_ID, int Strain_ID_R, int Strain_ID_D,  /* In use when defining the ODEs */
                                         double * Rate, Parameter_Table * Table );

void Preparing_Initial_System_Configuration( Parameter_Table * Table );  
void De_Allocating_Initial_System_Configuration ( Parameter_Table * Table );

void Time_Control_Upload ( Time_Control * Time, Parameter_Table * Table,
                                                          int I_Time, double Time_0, double Time_1 )
int Deterministic_Time_Dynamics( Parameter_Table * Table );
void Initial_Conditions_Numerical_Integration( Parameter_Table * Table, double * y_INI );
void Initial_Condition_Lower_Bound( Parameter_Table * Table, double * y_INI );							 
void Random_Initial_Condition( Parameter_Table * Table, double * y_INI );							 			
int numerical_Integration_Driver( Parameter_Table * Table, int j, double * Time_Current );
void JACOBIAN_Matrix(gsl_matrix *, const double *, double, Parameter_Table *);
double Per_Capita_Competition_Induced_Death_Calculation(Parameter_Table * Table, 
                                                        double * Y, int Strain_ID );
double Local_Population_Resources(int i, const double * Y, Parameter_Table * Table);		 
  


