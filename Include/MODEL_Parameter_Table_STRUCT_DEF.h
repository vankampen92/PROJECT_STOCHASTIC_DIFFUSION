/* Notice that Parameter_Table Struct depends and relies
   on the correct previous defintion of the structures
   Time_Control, Parameter_Model, and Parameter_Space
   but Parameter_Fitting type has not been defined yet. 
*/
typedef struct Parameter_Tableinfo
{
  /* The structure to fit parameters has not been defined yet. 
     I cannot do: 
  
     Parameter_Fitting * Fitting_Data; 
     
     I will get the same if I do: 
  */
  void * Fitting_Data;    /* When I need this to point to 
                             a Parameter_Fitting structure, 
                             I will have to do a cast!!! 
                          */
	
  int Focal_Resource;  /* Only for between-function communication purposes */

  /* Boolean Variable (see Time_Dependence_Control.c) */
  bool x_Bool; 

  /* Stochastic Realizations */
  int Realizations;

  /* * * Parameter Time Dependence Driven Dynamics * * */
  // #include "include.global.RainTemp.h" //or similar...
  // RainTemp * R_T;

  //"include.Trend_Control.global.h"
  Trend_Control * Tr;

  /* * * Time_Control * * * */
  Time_Control * T;

  /* * * Time_Dependence_Control * * * */
  Time_Dependence_Control * TDC;

  /* * * Parameter_Model * * * */
  #include "include.Parameter_Model.global.h"
  Parameter_Model * P;
    
  #include "include.Initial_Conditions.global.h"
  Parameter_Space * IC_Space;

  #include "include.Error_Control.global.h"
  Parameter_Space * E_Space; 

  /* * * Parameter_Space * * * */
  /* Total Number of Input Model Parameters (as in assign....c and boundary_*.c functions) */
  int MODEL_INPUT_PARAMETERS;
  char ** Name_Parameters; // Name_Parameters  : Name Model Input Parameters
  char ** Code_Parameters;
  char ** Symbol_Parameters;
  char ** Symbol_Parameters_Greek_CPGPLOT;
  
  double * Default_Vector_Parameters;

  /* What follows defines a parameter sub-space within
     the whole parameter space
  */
  int * Index;             // Vector defining the parameter subset to explore   */
  double * Vector_Parameters;
  //#include <include.PARAMETER_SPACE.global.h>
  Parameter_Space * S;
  /* * * * * * * * * */

  /* Total Number of Model Output Variables (in defintion_OutPut_Variables.c file) */
  int OUTPUT_VARIABLES_GENUINE; 
  int MODEL_OUTPUT_VARIABLES;
  char ** Output_Variable_Name;
  char ** Output_Variable_Symbol;
  double * Default_Vector_Output_Variables;
  /* Subset of Output Variables */
  /* What follows defines a subset of output variables that will be saved, represented, etc
   */
  int SUB_OUTPUT_VARIABLES;
  int * OUTPUT_VARIABLE_INDEX;
  double * Vector_Output_Variables;
  double ** Matrix_Output_Variables;
  /* * * * * * * * * */

  int SUM_LOCAL_STATE_VARIABLES; 
  int LOCAL_STATE_VARIABLES;
  int MODEL_STATE_VARIABLES;
  char ** Model_Variable_Name;
  char ** Model_Variable_Symbol;

  int K;
  int R;
  int A;
  int A_R; /* Index for Fertile Population */
  int A_A;
  int RA;
  int ARA;

  int * A_P;
  int * RA_P;
  int ** ARA_P; 

  double * Vector_Model_Variables; 
  double * Vector_Model_Variables_Time_0;

  double * Vector_Model_Variables_Stationarity;
  double ** Vector_Model_Variables_MultiStability;
  /*
     Evaluation of the largest (xmax) and the tiniest (xmin)
     numbers my machine can handle
  */
  float nr___x_min;
  float nr___x_MAX;

  double DETERMINISTIC_CASES; // Need by Error_Function.c 
  /* * * C P G   P L O T * * */

  int TOTAL_No_of_MODEL_PARAMETERS;
  
  double * Lambda_R;
  double * Delta_R; 
  
#if defined CPGPLOT_REPRESENTATION
#include <include.CPG.global.h>
  Parameter_CPGPLOT * CPG;
  Parameter_CPGPLOT * CPG_STO; 
#endif

  int * Vector_Model_Int_Variables;               /* Stochastic Dynamics */
  int * Vector_Model_Int_Variables_Time_0;        /* Stochastic Dynamics */

  double *** Metapop_Connectivity_Matrix; 
  Community ** Patch_System;                      /* Stochastic Dynamics */

  int TOTAL_No_of_EVENTS;                         /* Stochastic Dynamics */
  /* Number of common events that can occur to every Species: */
  int No_of_EVENTS;                               /* Stochastic Dynamics */

  /* TOTAL_No_of_EVENTS = No_of_EVENTS * No_of_RESOURCES if all species can undergo 
     exactly the same processes */
 
}Parameter_Table;

// void P_A_R_A_M_E_T_E_R___T_A_B_L_E___A_L_L_O_C( Parameter_Table * );

// void P_A_R_A_M_E_T_E_R___T_A_B_L_E___U_P_L_O_A_D( Parameter_Table *, int * );

// void Parameter_Values_into_Parameter_Table(Parameter_Table * );

// void values_TempDependency(double Temp, Parameter_Table *); ???

// void P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( Parameter_Table * );


