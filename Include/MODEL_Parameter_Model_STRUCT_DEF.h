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

#if defined CPGPLOT_REPRESENTATION
  Parameter_CPGPLOT * CPG;
#endif

  int K;

  int TOTAL_No_of_MODEL_PARAMETERS;

  double * Lambda_R;
  double * Delta_R; 
  
  int TOTAL_No_of_EVENTS;
  /* Total Number of Events within a patch, i.e., 25 * 4 = 100 */
  int No_of_EVENTS;

  double * Vector_Model_Variables; 
  double * Vector_Model_Variables_Time_0;
  int * Vector_Model_Int_Variables;               /* Stochastic Dynamics */
  int * Vector_Model_Int_Variables_Time_0;        /* Stochastic Dynamics */
  
  double *** Metapop_Connectivity_Matrix;
  
}Parameter_Model;
