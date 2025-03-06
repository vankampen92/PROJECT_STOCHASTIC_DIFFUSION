typedef struct Parameter_Modelinfo
{
  
  int MODEL_INPUT_PARAMETERS;
  int OUTPUT_VARIABLES_GENUINE; 
  int MODEL_OUTPUT_VARIABLES;
  int MODEL_STATE_VARIABLES;
  int LOCAL_STATE_VARIABLES;
  int SUM_LOCAL_STATE_VARIABLES; 
  
  /* * * Model Parameters  * * */
#include <include.Parameter_Model.global.h>

  /* * * Initial Conditions  * */
#include <include.Initial_Conditions.global.h>

  /* * * Error Model * * * */
#include "include.Error_Control.global.h"

  Trend_Control * Tr;

  Time_Control * Time; 

#if defined CPGPLOT_REPRESENTATION
  Parameter_CPGPLOT * CPG;
#endif

  int K;

  int TOTAL_No_of_MODEL_PARAMETERS;

  double * Lambda_R;
  double * Delta_R; 

  double * Nu_Consumers;   
  double * Theta_Consumers;

  double * y_R_i;
 
  double * Delta_RP;  /* Resource Propagule Death Rate (DIFFUSION_ECOEVO_PLANTS)              */

  double * Alpha_C;   /* Attack rates in Consumer-Resouce models                              */
                      /* Plasmid Reproduction Costs (DIFFUSION_ECO_PLASMIDS)                  */
  double * Nu_C;      /* Plasmid Resistance (DIFFUSION_ECO_PLASMIDS)                          */
                     
  double * Eta_RP;    /* Propagule Establishment Rate  (DIFFUSION_ECOEVO_PLANTS)              */
                      /* Strain-Specific Conjugation/Encounter Rates (DIFFUSION_ECO_PLASMIDS) */

  double * Beta_AP;   /* Propagule Production Rate of Adult Plants (DIFFUSION_ECOEVO_PLANTS)  */
                      /* Cell Division Rate of Single Bacterial Cell (DIFFUSION_ECO_PLASMIDS) */
  double * Delta_AP;  /* Death Rate of Adult Plants                (DIFFUSION_ECOEVO_PLANTS)  */
                      /* Per Capita Death Rates (DIFFUSION_ECO_PLASMIDS) */     
  
  double * Mu_RP;     /* Species-specific Diffusion/Movement Rate (DIFFUSION_ECOEVO_PLANTS)   */   
                      /* Bacterial Cell Diffusion/Movement Rates (DIFFUSION_ECO_PLASMIDS)     */
  double * Segregation_Error; /* Strain specific segregation errors (DIFFUSION_ECO_PLASMIDS)  */ 
  
  int TOTAL_No_of_EVENTS;
  /* Total Number of Events within a patch, i.e., 25 * 4 = 100 */
  int No_of_EVENTS;
  int TOTAL_GRAND_No_of_EVENTS;                   /* TOTAL_No_of_EVENTS * No_of_CELLS */ 

  double * Vector_Model_Variables; 
  double * Vector_Model_Variables_Time_0;
  int * Vector_Model_Int_Variables;               /* Stochastic Dynamics */
  int * Vector_Model_Int_Variables_Time_0;        /* Stochastic Dynamics */
  
  double *** Metapop_Connectivity_Matrix;

  treenode * Treeroot;                            /* Binary Sampling of a Discrete Distrubution */
  treenode ** Leaves;                             /* Leaves contain rates of individual events  */
  treenode *** Parent;                            /* Parent tree levels, but not the leaves     */ 

  int No_of_LEAVES;                               /* Imporant assert:                           */
  int No_of_TREE_LEVELS;                          /* (No_of_LEAVES == 2*No_of_TREE_LEVELS)      */

  /* Master Equation Numerical Integration */
  int TOTAL_No_of_RESOURCES;
  int TOTAL_No_of_CONSUMERS;
  int TOTAL_No_of_FREE_CONSUMERS_TIME_0;
  int TOTAL_No_of_HANDLING_CONSUMERS_TIME_0;
  int TOTAL_No_of_FREE_CONSUMERS;
  int TOTAL_No_of_HANDLING_CONSUMERS;
  
  void * Table; 
  
}Parameter_Model;
