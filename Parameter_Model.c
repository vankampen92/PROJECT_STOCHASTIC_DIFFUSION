#include "./Include/MODEL.h"

void Parameter_Model_Copy (Parameter_Model * P_Destination, Parameter_Model * P_Source)
{
  P_Destination->Mu_C       = P_Source->Mu_C;
  P_Destination->Mu         = P_Source->Mu;      /*  0 */

  P_Destination->Lambda_R_0 = P_Source->Lambda_R_0;
  P_Destination->Delta_R_0  = P_Source->Delta_R_0;
  P_Destination->Lambda_R_1 = P_Source->Lambda_R_1;
  P_Destination->Delta_R_1 = P_Source->Delta_R_1;
  P_Destination->K_R        = P_Source->K_R;

  P_Destination->Beta_R     = P_Source->Beta_R;        /* -H5 */

  P_Destination->Lambda_C_0 = P_Source->Lambda_C_0;     /* -H5 */
  P_Destination->Delta_C_0  = P_Source->Delta_C_0;      /* -H6 */
  P_Destination->Lambda_C_1 = P_Source->Lambda_C_1;     /* -H7 */
  P_Destination->Delta_C_1  = P_Source->Delta_C_1;      /* -H8 */

  P_Destination->Alpha_C_0  = P_Source->Alpha_C_0;      /* -H9 */
  P_Destination->Nu_C_0     = P_Source->Nu_C_0;         /* -H10 */

  P_Destination->Chi_C_0    = P_Source->Chi_C_0;        /* -H11 */
  P_Destination->Eta_C_0    = P_Source->Eta_C_0;        /* -H12 */

  P_Destination->N_E        = P_Source->N_E;   /* -H14 */ /* Number of Energy Levels */
  P_Destination->f          = P_Source->f;     /* -H15 */ /* Fecundity: No of Offspring Ind */ 
  P_Destination->i_0        = P_Source->i_0;   /* -H16 */ /* Energy Level at Maturity  */
  P_Destination->Beta_C     = P_Source->Beta_C;  /* -H17 */ /* Consummer Reproduction Rate */
  P_Destination->k_E        = P_Source->k_E;     /* -H18 */ /* 2*k_E resourse value energy units */
  P_Destination->Theta_C    = P_Source->Theta_C; /* -H19 */ /* Energy loss rate for maintenance */
  P_Destination->p_1        = P_Source->p_1;     /* -Hp1 */ /* 1st Cooperation probability  */ 
  P_Destination->p_2        = P_Source->p_2;    /* -Hp2 */  /* 2on Cooperation probability  */

  P_Destination->Eta_R      = P_Source->Eta_R;    /* -H20 */  /* Establishment rate  */ 

  P_Destination->No_of_IC = P_Source->No_of_IC;
  P_Destination->TYPE_of_INITIAL_CONDITION = P_Source->TYPE_of_INITIAL_CONDITION;
  P_Destination->INITIAL_TOTAL_POPULATION  = P_Source->INITIAL_TOTAL_POPULATION;
  P_Destination->RESCALING_INITIAL_TOTAL_POPULATION = P_Source->RESCALING_INITIAL_TOTAL_POPULATION;

  P_Destination->TYPE_of_ERROR_FUNCTION = P_Source->TYPE_of_ERROR_FUNCTION;
  P_Destination->No_of_ERROR_PARAMETERS = P_Source->No_of_ERROR_PARAMETERS;

  P_Destination->Error_Parameter_0 = P_Source->Error_Parameter_0;
  P_Destination->Error_Parameter_1 = P_Source->Error_Parameter_1;

  P_Destination->Err_0 = P_Source->Err_0;
  P_Destination->Err_1 = P_Source->Err_1;
  P_Destination->Err_2 = P_Source->Err_2;
  P_Destination->Err_3 = P_Source->Err_3;
  P_Destination->Err_4 = P_Source->Err_4;
  P_Destination->Err_5 = P_Source->Err_5;
  P_Destination->Err_6 = P_Source->Err_6;
  P_Destination->Err_7 = P_Source->Err_7;
  P_Destination->Err_8 = P_Source->Err_8;
  P_Destination->Err_9 = P_Source->Err_9;
  P_Destination->Err_10 = P_Source->Err_11;
  P_Destination->Err_11 = P_Source->Err_11;
  P_Destination->Err_12 = P_Source->Err_12;
  P_Destination->Err_13 = P_Source->Err_13;
  P_Destination->Err_14 = P_Source->Err_14;
  P_Destination->Err_15 = P_Source->Err_15;
  P_Destination->Err_16 = P_Source->Err_16;

  P_Destination->OUTPUT_VARIABLES_GENUINE = P_Source->OUTPUT_VARIABLES_GENUINE;/* Actual no of GENUINE (output) VARIABLES */
  P_Destination->MODEL_OUTPUT_VARIABLES = P_Source->MODEL_OUTPUT_VARIABLES;   /* Actual no of MODEL (output) VARIABLES */
  P_Destination->MODEL_INPUT_PARAMETERS = P_Source->MODEL_INPUT_PARAMETERS;  /* Actual no of MODEL (input) PARAMETERS */
  P_Destination->MODEL_STATE_VARIABLES  = P_Source->MODEL_STATE_VARIABLES;  /* Actual no of MODEL (state) VARIABLES  */
  P_Destination->LOCAL_STATE_VARIABLES  = P_Source->LOCAL_STATE_VARIABLES;  /* Actual no of MODEL (state) VARIABLES in every local population */
  P_Destination->SUM_LOCAL_STATE_VARIABLES  = P_Source->SUM_LOCAL_STATE_VARIABLES;  

  P_Destination->TYPE_of_MODEL        = P_Source->TYPE_of_MODEL;
  P_Destination->Growth_Function_Type = P_Source->Growth_Function_Type;

  //P_Destination->Time_Scale_Unit      = Table->T->Time_Scale_Unit;

  P_Destination->No_of_CELLS       = P_Source->No_of_CELLS;
  P_Destination->No_of_INDIVIDUALS = P_Source->No_of_INDIVIDUALS;
  P_Destination->No_of_CELLS_X     = P_Source->No_of_CELLS_X;
  P_Destination->No_of_CELLS_Y     = P_Source->No_of_CELLS_Y;
  P_Destination->TOTAL_GRAND_No_of_EVENTS = P_Source->TOTAL_GRAND_No_of_EVENTS;

  /* Ex: 8 (without acummulating variables) */
  P_Destination->No_of_NEIGHBORS            = P_Source->No_of_NEIGHBORS;
  P_Destination->TYPE_of_NETWORK            = P_Source->TYPE_of_NETWORK;
  P_Destination->TOTAL_No_of_EVENTS         = P_Source->TOTAL_No_of_EVENTS;
  /* Total Number of Events within a patch, i.e., 25 * 4 = 100 */
  P_Destination->No_of_EVENTS               = P_Source->No_of_EVENTS;
  /* Number of Events within an age class, i.e., 25           */

  P_Destination->Metapop_Connectivity_Matrix = P_Source->Metapop_Connectivity_Matrix;

  P_Destination->No_of_RESOURCES = P_Source->No_of_RESOURCES;
  P_Destination->No_of_PLASMIDS = P_Source->No_of_PLASMIDS;
  P_Destination->No_of_STRAINS  = P_Source->No_of_STRAINS;
  P_Destination->No_of_PROFILES = P_Source->No_of_PROFILES;

  /* Master Equation Numerical Integration */
  P_Destination->TOTAL_No_of_RESOURCES = P_Source->TOTAL_No_of_RESOURCES;
  P_Destination->TOTAL_No_of_CONSUMERS = P_Source->TOTAL_No_of_CONSUMERS;
  P_Destination->TOTAL_No_of_FREE_CONSUMERS_TIME_0 = P_Source->TOTAL_No_of_FREE_CONSUMERS_TIME_0;
  P_Destination->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0 = P_Source->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0;
  P_Destination->TOTAL_No_of_FREE_CONSUMERS = P_Source->TOTAL_No_of_FREE_CONSUMERS;
  P_Destination->TOTAL_No_of_HANDLING_CONSUMERS = P_Source->TOTAL_No_of_HANDLING_CONSUMERS;

  /* Pointers alineation between P and Table */
  P_Destination->Lambda_R = P_Source->Lambda_R;
  P_Destination->Delta_R  = P_Source->Delta_R; 

  P_Destination->Nu_Consumers    = P_Source->Nu_Consumers;
  P_Destination->Theta_Consumers = P_Source->Theta_Consumers;

  P_Destination->y_R_i     = P_Source->y_R_i;

  P_Destination->Delta_RP  = P_Source->Delta_RP;  /* Resource Propagule Death Rate (DIFFUSION_ECOEVO_PLANTS)              */

  P_Destination->Alpha_C  = P_Source->Alpha_C;   /* Attack rates in Consumer-Resouce models                              */
  P_Destination->Nu_C     = P_Source->Nu_C;      /* Plasmid Reproduction Costs (DIFFUSION_ECO_PLASMIDS) */

  P_Destination->Eta_RP   = P_Source->Eta_RP;                           /* Bacteria Conjugation Rate   */

  P_Destination->Beta_AP  = P_Source->Beta_AP;                          /* Bacteria Cell Division Rate */
  P_Destination->Delta_AP = P_Source->Delta_AP;                         /* Bacteria Death Rate         */

  P_Destination->Mu_RP    = P_Source->Mu_RP;                            /* Bacteria Diffusion Rates    */
  P_Destination->Segregation_Error = P_Source->Segregation_Error;                    
  /* End of pointer assignation */ 

  P_Destination->Time = P_Source->Time;
  P_Destination->Treeroot = P_Source->Treeroot; 
  P_Destination->Leaves = P_Source->Leaves; 
  P_Destination->Parent = P_Source->Parent; 
  P_Destination->No_of_LEAVES = P_Source->No_of_LEAVES;                                 
  P_Destination->No_of_TREE_LEVELS = P_Source->No_of_TREE_LEVELS;
  /* Imporant assert:                           */
  /* (No_of_LEAVES == 2*No_of_TREE_LEVELS)      */
}

void  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N ( Parameter_Table * Table,
							                                          Parameter_Model * P )
{
  /* This function transfer a subset of table parameters
     into the Parameter_Model structure. Parameter_Model parameters
     control the dynamical model. For back compatibilty, here 
      
     Model Parameter <--- Table Parameters, 

     when usually the first argument is the one being initialized (the destination) 
     and the second the one being copied (the source). 
  */

  P->Mu_C       = Table->Mu_C;
  P->Mu         = Table->Mu;       /*  0 */

  P->Lambda_R_0 = Table->Lambda_R_0;
  P->Delta_R_0  = Table->Delta_R_0;
  P->Lambda_R_1 = Table->Lambda_R_1;
  P->Delta_R_1  = Table->Delta_R_1;
  P->K_R        = Table->K_R;

  P->Beta_R     = Table->Beta_R;        /* -H5 */

  P->Lambda_C_0 = Table->Lambda_C_0;     /* -H5 */
  P->Delta_C_0  = Table->Delta_C_0;      /* -H6 */
  P->Lambda_C_1 = Table->Lambda_C_1;     /* -H7 */
  P->Delta_C_1  = Table->Delta_C_1;      /* -H8 */

  P->Alpha_C_0  = Table->Alpha_C_0;      /* -H9 */
  P->Nu_C_0     = Table->Nu_C_0;         /* -H10 */

  P->Chi_C_0    = Table->Chi_C_0;        /* -H11 */
  P->Eta_C_0    = Table->Eta_C_0;        /* -H12 */

  P->N_E        = Table->N_E;   /* -H14 */ /* Number of Energy Levels */
  P->f          = Table->f;     /* -H15 */ /* Fecundity: No of Offspring Ind */ 
  P->i_0        = Table->i_0;   /* -H16 */ /* Energy Level at Maturity  */
  P->Beta_C     = Table->Beta_C;  /* -H17 */ /* Consummer Reproduction Rate */
  P->k_E        = Table->k_E;     /* -H18 */ /* 2*k_E resourse value energy units */
  P->Theta_C    = Table->Theta_C; /* -H19 */ /* Energy loss rate for maintenance */
  P->p_1        = Table->p_1;     /* -Hp1 */ /* 1st Cooperation probability  */ 
  P->p_2        = Table->p_2;     /* -Hp2 */  /* 2on Cooperation probability  */ 

  P->Eta_R      = Table->Eta_R;    /* -H20 */  /* 2on Cooperation probability  */
  
  P->No_of_IC = Table->No_of_IC;
  P->TYPE_of_INITIAL_CONDITION = Table->TYPE_of_INITIAL_CONDITION;
  P->INITIAL_TOTAL_POPULATION  = Table->INITIAL_TOTAL_POPULATION;
  P->RESCALING_INITIAL_TOTAL_POPULATION = Table->RESCALING_INITIAL_TOTAL_POPULATION;

  P->TYPE_of_ERROR_FUNCTION = Table->TYPE_of_ERROR_FUNCTION;
  P->No_of_ERROR_PARAMETERS = Table->No_of_ERROR_PARAMETERS;
  P->Error_Parameter_0 = Table->Error_Parameter_0;
  P->Error_Parameter_1 = Table->Error_Parameter_1;

  P->Err_0 = Table->Err_0;
  P->Err_1 = Table->Err_1;
  P->Err_2 = Table->Err_2;
  P->Err_3 = Table->Err_3;
  P->Err_4 = Table->Err_4;
  P->Err_5 = Table->Err_5;
  P->Err_6 = Table->Err_6;
  P->Err_7 = Table->Err_7;
  P->Err_8 = Table->Err_8;
  P->Err_9 = Table->Err_9;
  P->Err_10 = Table->Err_11;
  P->Err_11 = Table->Err_11;
  P->Err_12 = Table->Err_12;
  P->Err_13 = Table->Err_13;
  P->Err_14 = Table->Err_14;
  P->Err_15 = Table->Err_15;
  P->Err_16 = Table->Err_16;

  P->OUTPUT_VARIABLES_GENUINE = Table->OUTPUT_VARIABLES_GENUINE;/* Actual no of MODEL (output) VARIABLES */
  P->MODEL_OUTPUT_VARIABLES = Table->MODEL_OUTPUT_VARIABLES;   /* Actual no of MODEL (output) VARIABLES */
  P->MODEL_INPUT_PARAMETERS = Table->MODEL_INPUT_PARAMETERS;  /* Actual no of MODEL (input) PARAMETERS */
  P->MODEL_STATE_VARIABLES  = Table->MODEL_STATE_VARIABLES;  /* Actual no of MODEL (state) VARIABLES  */
  P->LOCAL_STATE_VARIABLES  = Table->LOCAL_STATE_VARIABLES;  /* Actual no of MODEL (state) VARIABLES in every local population */
  P->SUM_LOCAL_STATE_VARIABLES  = Table->SUM_LOCAL_STATE_VARIABLES;  
  
  P->TYPE_of_MODEL        = Table->TYPE_of_MODEL;
  P->Growth_Function_Type = Table->Growth_Function_Type;

  //P->Time_Scale_Unit      = Table->T->Time_Scale_Unit;

  P->No_of_CELLS       = Table->No_of_CELLS;
  P->No_of_CELLS       = Table->No_of_CELLS;
  P->No_of_INDIVIDUALS = Table->No_of_INDIVIDUALS;
  P->No_of_CELLS_X     = Table->No_of_CELLS_X;
  P->No_of_CELLS_Y     = Table->No_of_CELLS_Y;
  P->TOTAL_GRAND_No_of_EVENTS = Table->TOTAL_GRAND_No_of_EVENTS;

  P->No_of_NEIGHBORS            = Table->No_of_NEIGHBORS;
  P->TYPE_of_NETWORK            = Table->TYPE_of_NETWORK;
  P->TOTAL_No_of_EVENTS         = Table->TOTAL_No_of_EVENTS;
  /* Total Number of Events within a patch, i.e., 18 * 4 = 72 */
  P->No_of_EVENTS               = Table->No_of_EVENTS;
  /* Number of Events within an age class, species or type , i.e., 18 */

  P->Metapop_Connectivity_Matrix = Table->Metapop_Connectivity_Matrix;

  P->No_of_RESOURCES = Table->No_of_RESOURCES;
  P->No_of_PLASMIDS  = Table->No_of_PLASMIDS;
  P->No_of_STRAINS   = Table->No_of_STRAINS;
  P->No_of_PROFILES  = Table->No_of_PROFILES;

  /* Master Equation Numerical Integration */
  P->TOTAL_No_of_RESOURCES = Table->TOTAL_No_of_RESOURCES;
  P->TOTAL_No_of_CONSUMERS = Table->TOTAL_No_of_CONSUMERS;
  P->TOTAL_No_of_FREE_CONSUMERS_TIME_0 = Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0;
  P->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0 = Table->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0;
  P->TOTAL_No_of_FREE_CONSUMERS = Table->TOTAL_No_of_FREE_CONSUMERS;
  P->TOTAL_No_of_HANDLING_CONSUMERS = Table->TOTAL_No_of_HANDLING_CONSUMERS;

  /* Pointers alineation between P and Table */
  P->Lambda_R = Table->Lambda_R;
  P->Delta_R  = Table->Delta_R; 

  P->Nu_Consumers    = Table->Nu_Consumers;
  P->Theta_Consumers = Table->Theta_Consumers;

  P->y_R_i     = Table->y_R_i;
  
  P->Delta_RP = Table->Delta_RP;      
  
  P->Alpha_C  = Table->Alpha_C;   /* Attack rates in Consumer-Resouce models                              */
  P->Nu_C     = Table->Nu_C;      /* Plasmid Reproduction Costs (DIFFUSION_ECO_PLASMIDS)                  */

  P->Eta_RP   = Table->Eta_RP;    /* Propagule Establishment Rate  (DIFFUSION_ECOEVO_PLANTS)              */
                                  /* Strain-Specific Conjugation/Encounter Rates (DIFFUSION_ECO_PLASMIDS) */ 
  P->Beta_AP  = Table->Beta_AP;                          /* Bacteria Cell Division Rate */
  P->Delta_AP = Table->Delta_AP;                         /* Bacteria Death Rate         */

  P->Mu_RP    = Table->Mu_RP;                            /* Bacteria Diffusion Rates    */
  P->Segregation_Error = Table->Segregation_Error;   
  /* End of pointer assignation */

  P->Time = Table->T;
  P->Treeroot = Table->Treeroot;
  P->Leaves = Table->Leaves;
  P->Parent = Table->Parent;  
  P->No_of_LEAVES = Table->No_of_LEAVES;                                 
  P->No_of_TREE_LEVELS = Table->No_of_TREE_LEVELS;
  /* Imporant assert:                           */
  /* (No_of_LEAVES == 2*No_of_TREE_LEVELS)      */

  P->Table = (void *)Table; 
}

void Parameter_Model_Copy_into_Parameter_Table (Parameter_Table * P_Destination, 
                                                Parameter_Model * P_Source)
{
  P_Destination->Mu_C       = P_Source->Mu_C;
  P_Destination->Mu         = P_Source->Mu;       /*  0 */

  P_Destination->Lambda_R_0 = P_Source->Lambda_R_0;
  P_Destination->Delta_R_0  = P_Source->Delta_R_0;
  P_Destination->Lambda_R_1 = P_Source->Lambda_R_1;
  P_Destination->Delta_R_1  = P_Source->Delta_R_1;
  P_Destination->K_R        = P_Source->K_R;

  P_Destination->Beta_R     = P_Source->Beta_R;        /* -H5 */

  P_Destination->Lambda_C_0 = P_Source->Lambda_C_0;     /* -H5 */
  P_Destination->Delta_C_0  = P_Source->Delta_C_0;      /* -H6 */
  P_Destination->Lambda_C_1 = P_Source->Lambda_C_1;     /* -H7 */
  P_Destination->Delta_C_1  = P_Source->Delta_C_1;      /* -H8 */

  P_Destination->Alpha_C_0  = P_Source->Alpha_C_0;      /* -H9 */
  P_Destination->Nu_C_0     = P_Source->Nu_C_0;         /* -H10 */

  P_Destination->Chi_C_0    = P_Source->Chi_C_0;        /* -H11 */
  P_Destination->Eta_C_0    = P_Source->Eta_C_0;        /* -H12 */

  P_Destination->N_E        = P_Source->N_E;   /* -H14 */ /* Number of Energy Levels */
  P_Destination->f          = P_Source->f;     /* -H15 */ /* Fecundity: No of Offspring Ind */ 
  P_Destination->i_0        = P_Source->i_0;   /* -H16 */ /* Energy Level at Maturity  */
  P_Destination->Beta_C     = P_Source->Beta_C;  /* -H17 */ /* Consummer Reproduction Rate */
  P_Destination->k_E        = P_Source->k_E;     /* -H18 */ /* 2*k_E resourse value energy units */
  P_Destination->Theta_C    = P_Source->Theta_C; /* -H19 */ /* Energy loss rate for maintenance */
  P_Destination->p_1        = P_Source->p_1;     /* -Hp1 */ /* 1st Cooperation probability  */ 
  P_Destination->p_2        = P_Source->p_2;    /* -Hp2 */  /* 2on Cooperation probability  */ 

  P_Destination->Eta_R      = P_Source->Eta_R;    /* -H20 */  /* 2on Cooperation probability  */ 
  
  P_Destination->No_of_IC = P_Source->No_of_IC;
  P_Destination->TYPE_of_INITIAL_CONDITION = P_Source->TYPE_of_INITIAL_CONDITION;
  P_Destination->INITIAL_TOTAL_POPULATION  = P_Source->INITIAL_TOTAL_POPULATION;
  P_Destination->RESCALING_INITIAL_TOTAL_POPULATION = P_Source->RESCALING_INITIAL_TOTAL_POPULATION;

  P_Destination->TYPE_of_ERROR_FUNCTION = P_Source->TYPE_of_ERROR_FUNCTION;
  P_Destination->No_of_ERROR_PARAMETERS = P_Source->No_of_ERROR_PARAMETERS;
  P_Destination->Error_Parameter_0 = P_Source->Error_Parameter_0;
  P_Destination->Error_Parameter_1 = P_Source->Error_Parameter_1;

  P_Destination->Err_0 = P_Source->Err_0;
  P_Destination->Err_1 = P_Source->Err_1;
  P_Destination->Err_2 = P_Source->Err_2;
  P_Destination->Err_3 = P_Source->Err_3;
  P_Destination->Err_4 = P_Source->Err_4;
  P_Destination->Err_5 = P_Source->Err_5;
  P_Destination->Err_6 = P_Source->Err_6;
  P_Destination->Err_7 = P_Source->Err_7;
  P_Destination->Err_8 = P_Source->Err_8;
  P_Destination->Err_9 = P_Source->Err_9;
  P_Destination->Err_10 = P_Source->Err_11;
  P_Destination->Err_11 = P_Source->Err_11;
  P_Destination->Err_12 = P_Source->Err_12;
  P_Destination->Err_13 = P_Source->Err_13;
  P_Destination->Err_14 = P_Source->Err_14;
  P_Destination->Err_15 = P_Source->Err_15;
  P_Destination->Err_16 = P_Source->Err_16;

  P_Destination->OUTPUT_VARIABLES_GENUINE = P_Source->OUTPUT_VARIABLES_GENUINE;/* Actual no of MODEL (output) VARIABLES */
  P_Destination->MODEL_OUTPUT_VARIABLES = P_Source->MODEL_OUTPUT_VARIABLES;   /* Actual no of MODEL (output) VARIABLES */
  P_Destination->MODEL_INPUT_PARAMETERS = P_Source->MODEL_INPUT_PARAMETERS;  /* Actual no of MODEL (input) PARAMETERS */
  P_Destination->MODEL_STATE_VARIABLES  = P_Source->MODEL_STATE_VARIABLES;  /* Actual no of MODEL (state) VARIABLES  */
  P_Destination->LOCAL_STATE_VARIABLES  = P_Source->LOCAL_STATE_VARIABLES;  /* Actual no of MODEL (state) VARIABLES in every local population */
  P_Destination->SUM_LOCAL_STATE_VARIABLES  = P_Source->SUM_LOCAL_STATE_VARIABLES;  

  P_Destination->TYPE_of_MODEL        = P_Source->TYPE_of_MODEL;
  P_Destination->Growth_Function_Type = P_Source->Growth_Function_Type;

  //P->Time_Scale_Unit      = Table->T->Time_Scale_Unit;

  P_Destination->No_of_CELLS       = P_Source->No_of_CELLS;
  P_Destination->No_of_INDIVIDUALS = P_Source->No_of_INDIVIDUALS;
  P_Destination->No_of_CELLS_X     = P_Source->No_of_CELLS_X;
  P_Destination->No_of_CELLS_Y     = P_Source->No_of_CELLS_Y;
  P_Destination->TOTAL_GRAND_No_of_EVENTS = P_Source->TOTAL_GRAND_No_of_EVENTS;

  /* Ex: 8 (without acummulating variables) */
  P_Destination->No_of_NEIGHBORS            = P_Source->No_of_NEIGHBORS;
  P_Destination->TYPE_of_NETWORK            = P_Source->TYPE_of_NETWORK;
  P_Destination->TOTAL_No_of_EVENTS         = P_Source->TOTAL_No_of_EVENTS;
  /* Total Number of Events within a patch, i.e., 18 * 4 = 72 */
  P_Destination->No_of_EVENTS               = P_Source->No_of_EVENTS;
  /* Number of Events within an age class, i.e., 18           */

  P_Destination->Metapop_Connectivity_Matrix = P_Source->Metapop_Connectivity_Matrix;

  P_Destination->No_of_RESOURCES = P_Source->No_of_RESOURCES;
  P_Destination->No_of_PLASMIDS = P_Source->No_of_PLASMIDS;
  P_Destination->No_of_STRAINS  = P_Source->No_of_STRAINS;
  P_Destination->No_of_PROFILES = P_Source->No_of_PROFILES;

  /* Master Equation Numerical Integration */
  P_Destination->TOTAL_No_of_RESOURCES = P_Source->TOTAL_No_of_RESOURCES;
  P_Destination->TOTAL_No_of_CONSUMERS = P_Source->TOTAL_No_of_CONSUMERS;
  P_Destination->TOTAL_No_of_FREE_CONSUMERS_TIME_0 = P_Source->TOTAL_No_of_FREE_CONSUMERS_TIME_0;
  P_Destination->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0 = P_Source->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0;
  P_Destination->TOTAL_No_of_FREE_CONSUMERS = P_Source->TOTAL_No_of_FREE_CONSUMERS;
  P_Destination->TOTAL_No_of_HANDLING_CONSUMERS = P_Source->TOTAL_No_of_HANDLING_CONSUMERS;

  /* Pointers alineation between P and Table */
  P_Destination->Lambda_R = P_Source->Lambda_R;
  P_Destination->Delta_R  = P_Source->Delta_R; 

  P_Destination->Nu_Consumers    = P_Source->Nu_Consumers;
  P_Destination->Theta_Consumers = P_Source->Theta_Consumers;

  P_Destination->y_R_i     = P_Source->y_R_i;

  P_Destination->Delta_RP  = P_Source->Delta_RP;  /* Resource Propagule Death Rate (DIFFUSION_ECOEVO_PLANTS)              */

  P_Destination->Alpha_C  = P_Source->Alpha_C;   /* Attack rates in Consumer-Resouce models                              */
  P_Destination->Nu_C     = P_Source->Nu_C;      /* Plasmid Reproduction Costs (DIFFUSION_ECO_PLASMIDS) */

  P_Destination->Eta_RP   = P_Source->Eta_RP;                           /* Bacteria Conjugation Rate   */

  P_Destination->Beta_AP  = P_Source->Beta_AP;                          /* Bacteria Cell Division Rate */
  P_Destination->Delta_AP = P_Source->Delta_AP;                         /* Bacteria Death Rate         */

  P_Destination->Mu_RP    = P_Source->Mu_RP;                            /* Bacteria Diffusion Rates    */
  P_Destination->Segregation_Error = P_Source->Segregation_Error;                    
  /* End of pointer assignation */
 
  P_Destination->T = P_Source->Time;
  P_Destination->Treeroot = P_Source->Treeroot; 
  P_Destination->Leaves = P_Source->Leaves; 
  P_Destination->Parent = P_Source->Parent; 
  P_Destination->No_of_LEAVES = P_Source->No_of_LEAVES;                                 
  P_Destination->No_of_TREE_LEVELS = P_Source->No_of_TREE_LEVELS;
  /* Imporant assert:                           */
  /* (No_of_LEAVES == 2*No_of_TREE_LEVELS)      */
}

void Parameter_Table_Copy_into_Parameter_Model (Parameter_Model * P_Destination, 
                                                Parameter_Table * P_Source)
{
  P_Destination->Mu_C       = P_Source->Mu_C;
  P_Destination->Mu         = P_Source->Mu;       /*  0 */

  P_Destination->Lambda_R_0 = P_Source->Lambda_R_0;
  P_Destination->Delta_R_0  = P_Source->Delta_R_0;
  P_Destination->Lambda_R_1 = P_Source->Lambda_R_1;
  P_Destination->Delta_R_1  = P_Source->Delta_R_1;
  P_Destination->K_R        = P_Source->K_R;

  P_Destination->Beta_R     = P_Source->Beta_R;        /* -H5 */

  P_Destination->Lambda_C_0 = P_Source->Lambda_C_0;     /* -H5 */
  P_Destination->Delta_C_0  = P_Source->Delta_C_0;      /* -H6 */
  P_Destination->Lambda_C_1 = P_Source->Lambda_C_1;     /* -H7 */
  P_Destination->Delta_C_1  = P_Source->Delta_C_1;      /* -H8 */

  P_Destination->Alpha_C_0  = P_Source->Alpha_C_0;      /* -H9 */
  P_Destination->Nu_C_0     = P_Source->Nu_C_0;         /* -H10 */

  P_Destination->Chi_C_0    = P_Source->Chi_C_0;        /* -H11 */
  P_Destination->Eta_C_0    = P_Source->Eta_C_0;        /* -H12 */

  P_Destination->N_E        = P_Source->N_E;   /* -H14 */ /* Number of Energy Levels */
  P_Destination->f          = P_Source->f;     /* -H15 */ /* Fecundity: No of Offspring Ind */ 
  P_Destination->i_0        = P_Source->i_0;   /* -H16 */ /* Energy Level at Maturity  */
  P_Destination->Beta_C     = P_Source->Beta_C;  /* -H17 */ /* Consummer Reproduction Rate */
  P_Destination->k_E        = P_Source->k_E;     /* -H18 */ /* 2*k_E resourse value energy units */
  P_Destination->Theta_C    = P_Source->Theta_C; /* -H19 */ /* Energy loss rate for maintenance */
  P_Destination->p_1        = P_Source->p_1;     /* -Hp1 */ /* 1st Cooperation probability  */ 
  P_Destination->p_2        = P_Source->p_2;    /* -Hp2 */  /* 2on Cooperation probability  */ 

  P_Destination->Eta_R      = P_Source->Eta_R;    /* -H20 */  /* 2on Cooperation probability  */ 
  
  P_Destination->No_of_IC = P_Source->No_of_IC;
  P_Destination->TYPE_of_INITIAL_CONDITION = P_Source->TYPE_of_INITIAL_CONDITION;
  P_Destination->INITIAL_TOTAL_POPULATION  = P_Source->INITIAL_TOTAL_POPULATION;
  P_Destination->RESCALING_INITIAL_TOTAL_POPULATION = P_Source->RESCALING_INITIAL_TOTAL_POPULATION;

  P_Destination->TYPE_of_ERROR_FUNCTION = P_Source->TYPE_of_ERROR_FUNCTION;
  P_Destination->No_of_ERROR_PARAMETERS = P_Source->No_of_ERROR_PARAMETERS;
  P_Destination->Error_Parameter_0 = P_Source->Error_Parameter_0;
  P_Destination->Error_Parameter_1 = P_Source->Error_Parameter_1;

  P_Destination->Err_0 = P_Source->Err_0;
  P_Destination->Err_1 = P_Source->Err_1;
  P_Destination->Err_2 = P_Source->Err_2;
  P_Destination->Err_3 = P_Source->Err_3;
  P_Destination->Err_4 = P_Source->Err_4;
  P_Destination->Err_5 = P_Source->Err_5;
  P_Destination->Err_6 = P_Source->Err_6;
  P_Destination->Err_7 = P_Source->Err_7;
  P_Destination->Err_8 = P_Source->Err_8;
  P_Destination->Err_9 = P_Source->Err_9;
  P_Destination->Err_10 = P_Source->Err_11;
  P_Destination->Err_11 = P_Source->Err_11;
  P_Destination->Err_12 = P_Source->Err_12;
  P_Destination->Err_13 = P_Source->Err_13;
  P_Destination->Err_14 = P_Source->Err_14;
  P_Destination->Err_15 = P_Source->Err_15;
  P_Destination->Err_16 = P_Source->Err_16;

  P_Destination->OUTPUT_VARIABLES_GENUINE = P_Source->OUTPUT_VARIABLES_GENUINE;/* Actual no of MODEL (output) VARIABLES */
  P_Destination->MODEL_OUTPUT_VARIABLES = P_Source->MODEL_OUTPUT_VARIABLES;   /* Actual no of MODEL (output) VARIABLES */
  P_Destination->MODEL_INPUT_PARAMETERS = P_Source->MODEL_INPUT_PARAMETERS;  /* Actual no of MODEL (input) PARAMETERS */
  P_Destination->MODEL_STATE_VARIABLES  = P_Source->MODEL_STATE_VARIABLES;  /* Actual no of MODEL (state) VARIABLES  */
  P_Destination->LOCAL_STATE_VARIABLES  = P_Source->LOCAL_STATE_VARIABLES;  /* Actual no of MODEL (state) VARIABLES in every local population */
  P_Destination->SUM_LOCAL_STATE_VARIABLES  = P_Source->SUM_LOCAL_STATE_VARIABLES;  

  P_Destination->TYPE_of_MODEL        = P_Source->TYPE_of_MODEL;
  P_Destination->Growth_Function_Type = P_Source->Growth_Function_Type;
  P_Destination->TYPE_of_NETWORK      = P_Source->TYPE_of_NETWORK;

  //P->Time_Scale_Unit      = Table->T->Time_Scale_Unit;

  P_Destination->No_of_CELLS       = P_Source->No_of_CELLS;
  P_Destination->No_of_INDIVIDUALS = P_Source->No_of_INDIVIDUALS;
  P_Destination->No_of_CELLS_X     = P_Source->No_of_CELLS_X;
  P_Destination->No_of_CELLS_Y     = P_Source->No_of_CELLS_Y;
  P_Destination->No_of_NEIGHBORS            = P_Source->No_of_NEIGHBORS;
  P_Destination->TOTAL_No_of_EVENTS         = P_Source->TOTAL_No_of_EVENTS;
  P_Destination->TOTAL_GRAND_No_of_EVENTS = P_Source->TOTAL_GRAND_No_of_EVENTS;

  /* Total Number of Events within a patch, i.e., 18 * 4 = 72 */
  P_Destination->No_of_EVENTS               = P_Source->No_of_EVENTS;
  /* Number of Events within an age class, i.e., 18           */

  P_Destination->Metapop_Connectivity_Matrix = P_Source->Metapop_Connectivity_Matrix;

  P_Destination->No_of_RESOURCES = P_Source->No_of_RESOURCES;
  P_Destination->No_of_PLASMIDS = P_Source->No_of_PLASMIDS;
  P_Destination->No_of_STRAINS  = P_Source->No_of_STRAINS;
  P_Destination->No_of_PROFILES = P_Source->No_of_PROFILES;
  
  /* Master Equation Numerical Integration */
  P_Destination->TOTAL_No_of_RESOURCES = P_Source->TOTAL_No_of_RESOURCES;
  P_Destination->TOTAL_No_of_CONSUMERS = P_Source->TOTAL_No_of_CONSUMERS;
  P_Destination->TOTAL_No_of_FREE_CONSUMERS_TIME_0 = P_Source->TOTAL_No_of_FREE_CONSUMERS_TIME_0;
  P_Destination->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0 = P_Source->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0;
  P_Destination->TOTAL_No_of_FREE_CONSUMERS = P_Source->TOTAL_No_of_FREE_CONSUMERS;
  P_Destination->TOTAL_No_of_HANDLING_CONSUMERS = P_Source->TOTAL_No_of_HANDLING_CONSUMERS;

  /* Pointers alineation between P and Table */
  P_Destination->Lambda_R = P_Source->Lambda_R;
  P_Destination->Delta_R  = P_Source->Delta_R; 

  P_Destination->Nu_Consumers    = P_Source->Nu_Consumers;
  P_Destination->Theta_Consumers = P_Source->Theta_Consumers;

  P_Destination->y_R_i     = P_Source->y_R_i;

  P_Destination->Delta_RP  = P_Source->Delta_RP;  /* Resource Propagule Death Rate (DIFFUSION_ECOEVO_PLANTS)              */

  P_Destination->Alpha_C  = P_Source->Alpha_C;   /* Attack rates in Consumer-Resouce models                              */
  P_Destination->Nu_C     = P_Source->Nu_C;      /* Plasmid Reproduction Costs (DIFFUSION_ECO_PLASMIDS) */

  P_Destination->Eta_RP   = P_Source->Eta_RP;                           /* Bacteria Conjugation Rate   */

  P_Destination->Beta_AP  = P_Source->Beta_AP;                          /* Bacteria Cell Division Rate */
  P_Destination->Delta_AP = P_Source->Delta_AP;                         /* Bacteria Death Rate         */

  P_Destination->Mu_RP    = P_Source->Mu_RP;                            /* Bacteria Diffusion Rates    */
  P_Destination->Segregation_Error = P_Source->Segregation_Error;                    
  /* End of pointer assignation */
 

  P_Destination->Time = P_Source->T;
  P_Destination->Treeroot = P_Source->Treeroot;  
  P_Destination->Leaves = P_Source->Leaves; 
  P_Destination->Parent = P_Source->Parent;
  P_Destination->No_of_LEAVES = P_Source->No_of_LEAVES;                                 
  P_Destination->No_of_TREE_LEVELS = P_Source->No_of_TREE_LEVELS;
  /* Imporant assert:                           */
  /* (No_of_LEAVES == 2*No_of_TREE_LEVELS)      */
}

void Vector_Entries_into_Parameter_Model ( const gsl_vector * X, Parameter_Model * P,
					                                  int * Parameter_Index, int No_of_PARAMETERS )
{
  int i;
  int key;
  double value;

  for( i=0; i<No_of_PARAMETERS; i++) {
    key = Parameter_Index[i];
    value = gsl_vector_get(X, i);
    Vector_Entry_into_Parameter_Model ( value, key, P );
  }
}

void Vector_Entry_into_Parameter_Model ( double value, int key, Parameter_Model * P )
{

  switch (key) {

  case  0: P->Mu = value;
    break;
  case  1: P->No_of_INDIVIDUALS = (int)value;
    break;
  case  2: P->No_of_CELLS       = (int)value;
    break;
  case  3: P->No_of_CELLS_X     = (int)value;
    break;
  case  4: P->No_of_CELLS_Y     = (int)value;
    break;
  case  5: P->No_of_RESOURCES     = (int)value;
    break;

  case  6: P->Lambda_R_0 = value;
    break;
  case  7: P->Delta_R_0  = value;
    break;
  case  8: P->Lambda_R_1 = value;
    break;
  case  9: P->Delta_R_1  = value;
    break;
  case 10: P->K_R        = (int)value;              /* Resource Carrying Capacity */
    break;
  case 11: P->Beta_R = value;        /* -H4 */
    break;

  case 12: P->Lambda_C_0 = value;     /* -H5 */
    break;
  case 13: P->Delta_C_0 = value;      /* -H6 */
    break;
  case 14: P->Lambda_C_1 = value;     /* -H7 */
    break;
  case 15: P->Delta_C_1 = value;      /* -H8 */
    break;

  case 16: P->Alpha_C_0 = value;      /* -H9 */
    break;
  case 17: P->Nu_C_0 = value;         /* -H10 */
    break;

  case 18: P->Chi_C_0 = value;        /* -H11 */
    break;
  case 19: P->Eta_C_0 = value;        /* -H12 */
    break;

  case 20: P->Mu_C = value;           /* -H13 */
    break;

  case 21: P->N_E        = value; //N_E;  /* -H14 */ /* Number of Energy Levels */
    break;
    
  case 22: P->f          = value; //f;    /* -H15 */ /* Fecundity: Number of Offspring  */
    break;
    
  case 23: P->i_0        = value; //i_0;  /* -H16 */ /* Energy Level at Maturity  */
    break;
    
  case 24: P->Beta_C     = value; //Beta_C;  /* -H17 */ /* Consummer Reproduction Rate */
    break;
    
  case 25: P->k_E        = value; /* -H18 */ /* 2* k_E is the resourse value in energy units */
    break;
    
  case 26: P->Theta_C    = value; /* -H19 */ /* Energy loss rate for maintenance */
    break;
    
  case 27: P->p_1        = value; /* Cooperation probability 1st position in the triplet */
    break;
    
  case 28: P->p_2        = value; /* Cooperation probability 2on position in the triplet */  
    break;

  case 29: P->Eta_R     = value; /* Propagule Establishment Rate  */  
    break;

  default:
    printf(".... INVALID PARAMETER KEY (key = %d)\n", key);
    printf(" The maximum number of parameters is Number_PAR_MAX\n");
    printf(" The permited number of keys go from 0, to %d\n", MODEL_PARAMETERS_MAXIMUM-1);

    exit(0);
  }
}

void Parameter_Model_into_Vector_Entries ( Parameter_Model * P, gsl_vector * X,
					   int * Parameter_Index, int No_of_PARAMETERS )
{
  int i;
  int key;
  double value;

  for( i=0; i<No_of_PARAMETERS; i++) {
    key = Parameter_Index[i];
    value = Parameter_Model_into_Vector_Entry( key, P );
    gsl_vector_set(X, i, value);
  }
}

double Parameter_Model_into_Vector_Entry ( int key, Parameter_Model * P )
{
  double value;

  switch (key) {

  case  0: value = P->Mu;
    break;
  case  1: value = (double)P->No_of_INDIVIDUALS;
    break;
  case  2: value = (double)P->No_of_CELLS;
    break;
  case  3: value = (double)P->No_of_CELLS_X;
    break;
  case  4: value = (double)P->No_of_CELLS_Y;
    break;
  case  5: value = (double)P->No_of_RESOURCES;
    break;

  case  6: value = P->Lambda_R_0;
      break;
  case  7: value = P->Delta_R_0;
    break;
  case  8: value = P->Lambda_R_1;
    break;
  case  9: value = P->Delta_R_1;
    break;
  case 10: value = (double)P->K_R;   /* Resource Carrying Capacity */
    break;
  case 11: value = P->Beta_R;        /* -H4 */
    break;

  case 12: value = P->Lambda_C_0;     /* -H5 */
    break;
  case 13: value = P->Delta_C_0;      /* -H6 */
    break;
  case 14: value = P->Lambda_C_1;     /* -H7 */
    break;
  case 15: value = P->Delta_C_1;      /* -H8 */
    break;

  case 16: value = P->Alpha_C_0;      /* -H9 */
    break;
  case 17: value = P->Nu_C_0;         /* -H10 */
    break;

  case 18: value = P->Chi_C_0;        /* -H11 */
    break;
  case 19: value = P->Eta_C_0;        /* -H12 */
    break;

  case 20: value = P->Mu_C;           /* -H13 */
    break;

  case 21: value = P->N_E        ; //N_E;  /* -H14 */ /* Number of Energy Levels */
    break;
    
  case 22: value = P->f          ; //f;    /* -H15 */ /* Fecundity: Number of Offspring  */
    break;
    
  case 23: value = P->i_0        ; //i_0;  /* -H16 */ /* Energy Level at Maturity  */
    break;
    
  case 24: value = P->Beta_C     ; //Beta_C;  /* -H17 */ /* Consummer Reproduction Rate */
    break;
    
  case 25: value = P->k_E        ; /* -H18 */ /* 2* k_E is the resourse value in energy units */
    break;
    
  case 26: value = P->Theta_C    ; /* -H19 */ /* Energy loss rate for maintenance */
    break;
    
  case 27: value = P->p_1        ; /* Cooperation probability 1st position in the triplet */
    break;
    
  case 28: value = P->p_2        ; /* Cooperation probability 2on position in the triplet */  
    break;

  case 29: value = P->Eta_R      ; /* Propagule Establishment Rate  */  
    break;

  default:
      printf(".... INVALID PARAMETER KEY (key = %d)\n", key);
      printf(" The maximum number of parameters is Number_PAR_MAX\n");
      printf(" The permited number of keys go from 0, to %d\n", MODEL_PARAMETERS_MAXIMUM-1);

      exit(0);
  }

  return(value);
}

void P_A_R_A_M_E_T_E_R___M_O_D_E_L___F_R_E_E( Parameter_Model * P)
{
  free ( P );
}
