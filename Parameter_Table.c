#include <MODEL.h>

#if defined CPGPLOT_REPRESENTATION
#include <include.CPG.extern.h>
#endif

#include <include.Parameter_Model.extern.h>
#include <include.Output_Variables.extern.h>
#include <include.Parameter_Space.extern.h>
#include <include.Time_Control.extern.h>
#include <include.Initial_Conditions.extern.h>
#include <include.Error_Control.extern.h>
#include <include.Type_of_Model.extern.h>

// #include "./Include/include.Stochastic_Control.extern.h"
extern int Realizations;

void P_A_R_A_M_E_T_E_R___T_A_B_L_E___A_L_L_O_C( Parameter_Table * Table )
{
  int i, j, a;
  int MODEL_OUTPUT_VARIABLES;   /* Actual no of MODEL (output) VARIABLES */
  int MODEL_INPUT_PARAMETERS;  /* Actual no of MODEL (input) PARAMETERS */

  MODEL_INPUT_PARAMETERS = MODEL_PARAMETERS_MAXIMUM;
  MODEL_OUTPUT_VARIABLES = OUTPUT_VARIABLES_MAXIMUM;
  // int MODEL_STATE_VARIABLES;  /* Actual no of MODEL (state) VARIABLES 
  /* Model parameters are all input parameters defining and controling
     model dynamics (both through simulation and mathematically, in case
     the model has a mathematical description such as a system of ODEs)
  */
  Table->Code_Parameters = (char **)malloc( MODEL_INPUT_PARAMETERS * sizeof(char *) );
  Table->Symbol_Parameters = (char **)malloc( MODEL_INPUT_PARAMETERS * sizeof(char *) );
  Table->Symbol_Parameters_Greek_CPGPLOT = (char **)malloc( MODEL_INPUT_PARAMETERS * sizeof(char *));
  Table->Name_Parameters = (char **)malloc( MODEL_INPUT_PARAMETERS * sizeof(char *) );
  for (i=0; i<MODEL_INPUT_PARAMETERS; i++){
    Table->Code_Parameters[i] = (char *)malloc( 100 * sizeof(char) );
    Table->Name_Parameters[i] = (char *)malloc( 100 * sizeof(char) );
    Table->Symbol_Parameters[i] = (char *)malloc( 100 * sizeof(char) );
    Table->Symbol_Parameters_Greek_CPGPLOT[i] = (char *)malloc( 100 * sizeof(char) );
  }
  Table->Default_Vector_Parameters = (double *)malloc( MODEL_INPUT_PARAMETERS * sizeof(double) );

  Table->Index = (int *)malloc( MODEL_PARAMETERS_MAXIMUM * sizeof(int) );
  Table->Vector_Parameters = (double *)malloc( MODEL_PARAMETERS_MAXIMUM * sizeof(double) );

  /* Output Variables are any measure of the state model variables of any function of
     these at any given time:
     MODEL_OUTPUT_VARIABLES is the Total Number of Total Output Variables
     (as it appears in defintion_OutPut_Variables.c file)( )
   */
  Table->Output_Variable_Name = (char **)malloc( MODEL_OUTPUT_VARIABLES * sizeof(char *) );
  Table->Output_Variable_Symbol = (char **)malloc( MODEL_OUTPUT_VARIABLES * sizeof(char *) );
  for (i=0; i<MODEL_OUTPUT_VARIABLES; i++){
    Table->Output_Variable_Name[i] = (char *)malloc( 100 * sizeof(char) );
    Table->Output_Variable_Symbol[i] = (char *)malloc( 100 * sizeof(char) );
  }
  Table->OUTPUT_VARIABLE_INDEX = (int *)malloc( SUB_OUTPUT_VARIABLES * sizeof(int) );
  Table->Vector_Output_Variables = (double *)malloc( SUB_OUTPUT_VARIABLES * sizeof(double) );
  Table->Matrix_Output_Variables = (double **)malloc( SUB_OUTPUT_VARIABLES * sizeof(double *) );
  for (i=0; i<SUB_OUTPUT_VARIABLES; i++)
    Table->Matrix_Output_Variables[i] = (double *)malloc( I_Time * sizeof(double) );

  Table->Default_Vector_Output_Variables = (double *)malloc( MODEL_OUTPUT_VARIABLES * sizeof(double) );

  Table->Vector_Model_Variables_MultiStability =(double **)calloc( 3, sizeof(double *) );

  /* Model Variables represent only the state model variables, i.e., the set of variables
     completely defining the state of the system at any given time. Model Variables can be,
     of course, also output variables. They can be defined as output variables
     whose definition can be found at Definition_Output_Variables.c */

  /* Some implementations of this code require to alloc memmory according to
     a number of state variables that can change dynamically. In that case, this
     6 lines of code should be commented out. Notice alse the corresponding lines
     in void P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( Parameter_Table * T )

     For most calculations, I opted for allocating this part of the data structure
     only if needed (see ./MODEL_CALCULATIONS/TIME_EVO_DETERMINISTIC/MODEL.c) regardless
     the number of MODEL STATE VARIABLES changing or not during execution. Essentially,
     the numerical integration of a system of ODEs relies on Table->Vector_Model_Variables.
     Therefore, a driver, usaully called void MODEL(Parameter_Table *), is in charge
     of doing this allocation previous to calling the numerical integration function,
     and deallocating this memmory right after.
  */
  // Table->Model_Variable_Name = (char **)malloc( MODEL_STATE_VARIABLES * sizeof(char *) );
  // for (i=0; i<MODEL_STATE_VARIABLES; i++){
  //   Table->Model_Variable_Name[i] = (char *)malloc( 100 * sizeof(char) );
  // }
  // Table->Vector_Model_Variables = (double *)malloc( MODEL_STATE_VARIABLES * sizeof(double) );
  // Table->Vector_Model_Variables_Stationarity = (double *)malloc( MODEL_STATE_VARIABLES * sizeof(double) );

  if(     TYPE_of_NETWORK == 0) No_of_NEIGHBORS = No_of_CELLS - 1;
  else if(TYPE_of_NETWORK == 1) No_of_NEIGHBORS = 4;               /* Von Neumann */
  else {
    printf(" No more network types have been defined so far\n");
    printf(" The program will exit\n");
    exit(0); 
  }
  
  /* BEGIN: Allocating and Setting up Connectivity Matrix */
  Table->Metapop_Connectivity_Matrix = (double ***)calloc(No_of_SPECIES,
							  sizeof(double **) );
  for(a=0; a<No_of_SPECIES; a++) { 
    Table->Metapop_Connectivity_Matrix[a] = (double **)calloc(No_of_CELLS,
							      sizeof(double *) );
    for(i=0; i<No_of_CELLS; i++) {
      Table->Metapop_Connectivity_Matrix[a][i] = (double *)calloc(No_of_NEIGHBORS+1,
								  sizeof(double) );
    }
  }
  /* END: Allocating and Setting up Connectivity Matrix */
}

void P_A_R_A_M_E_T_E_R___T_A_B_L_E___U_P_L_O_A_D( Parameter_Table * Table, int * Index_Output_Variables )
{
  int i, j, a;
  
  /* Stochastic Realizations */
  Table->Realizations = Realizations;

  /* Total number of potential input paramters */
  Table->MODEL_INPUT_PARAMETERS = MODEL_PARAMETERS_MAXIMUM;

  Table->No_of_CELLS       = No_of_CELLS; 
  Table->No_of_INDIVIDUALS = No_of_INDIVIDUALS;
  Table->No_of_CELLS_X     = No_of_CELLS_X;
  Table->No_of_CELLS_Y     = No_of_CELLS_Y;
  
  /* Parameter Model Upload */
  Parameter_Values_into_Parameter_Table(Table);

  /* Type of Model upload  */
  Table->TYPE_of_MODEL = TYPE_of_MODEL;  assert_right_model_definition( Table );
  Model_Variables_Code_into_Parameter_Table (Table);
  /* Total number of potential state variables */
  /* Total number of potential state variables */
  Table->MODEL_STATE_VARIABLES = Table->K + 1;
  /* Total number of potential output variables */
  Table->MODEL_OUTPUT_VARIABLES = OUTPUT_VARIABLES_GENUINE + Table->MODEL_STATE_VARIABLES;

  /* Total number of actual model output variables */
  Table->SUB_OUTPUT_VARIABLES   = SUB_OUTPUT_VARIABLES;
  /* Parameter Space Upload: PARAMETER SPACE             */
  /* Total number of actual input paramters              */
  /* Table->No_of_PARAMETERS = MODEL_PARAMETERS_MAXIMUM; */
  /* Table->A_n = A_n;                                   */
  /* Table->A_d = A_d;                                   */
  /* Table->No_of_POINTS    = No_of_POINTS;              */
  /* Table->Input_Parameter = Input_Parameter;           */
  /* Table->Value_0         = Value_0;                   */
  /* Table->Value_1         = Value_1;                   */

  /// Setting MODEL INPUT PARAMETERS up !!!
  Table->Growth_Function_Type = Growth_Function_Type;
  /* BEGIN: Parameter default values into vector structure */
  for(i = 0; i<Table->MODEL_INPUT_PARAMETERS; i++){
    Table->Default_Vector_Parameters[i] = AssignStructValue_to_VectorEntry(i, Table);
  }
  /*   END: Parameter default values into vector structure */

  /* BEGIN: Names and codes for model parameters */
  for(i = 0; i<Table->MODEL_INPUT_PARAMETERS; i++){
    AssignLabel_to_Model_Parameters(i, Table->Name_Parameters[i], Table);
    AssignCodes_to_Model_Parameters(i, Table->Code_Parameters[i], Table);
    AssignSymbol_to_Model_Parameters(i, Table->Symbol_Parameters[i], Table);
    AssignCPGPLOT_Symbol_to_Model_Parameters(i, Table->Symbol_Parameters_Greek_CPGPLOT[i], Table); 
  }
  /*   END: Names and codes for model parameter  */

  /* Default model input parameters */
  for (i=0; i < MODEL_PARAMETERS_MAXIMUM; i++) Table->Index[i] = i;
  
  // This assignation will be overwritten when Parameter_Space structure is setup

  /// Setting MODEL STATE VARIABLES up !!!  
  int MODEL_STATE_VARIABLES    = Table->K + 1;
  Table->Model_Variable_Name   = (char **)malloc( MODEL_STATE_VARIABLES * sizeof(char *) );
  Table->Model_Variable_Symbol = (char **)malloc( MODEL_STATE_VARIABLES * sizeof(char *) );
  for (i=0; i<MODEL_STATE_VARIABLES; i++){
     Table->Model_Variable_Name[i] = (char *)malloc( 100 * sizeof(char) );
     Table->Model_Variable_Symbol[i] = (char *)malloc( 100 * sizeof(char) );
  }
  /* BEGIN: Names and symbols for model state variables */
  for(i = 0; i<Table->MODEL_STATE_VARIABLES; i++){
    AssignLabel_to_Model_Variables(i, Table->Model_Variable_Name[i], Table);
    AssignSymbol_to_Model_Variables(i, Table->Model_Variable_Symbol[i], Table);
  }
  /*   END: Names and codes for model variables  */

  /// Setting MODEL OUTPUT VARIABLES up !!!  
  /* BEGIN: Names and symbol for output variables           */
  for(i = 0; i<Table->MODEL_OUTPUT_VARIABLES; i++){
    AssignLabel_to_Output_Variables(i, Table->Output_Variable_Name[i], Table);
    // AssignSymbol_to_Output_Variables(i, Table->Output_Variable_Symbol[i], Table);
    AssignCPGPLOT_Symbol_to_Output_Variables(i, Table->Output_Variable_Symbol[i], Table);
  }
  /*   END: Names for output variables */
  for(i=0; i < SUB_OUTPUT_VARIABLES; i++)
    Table->OUTPUT_VARIABLE_INDEX[i] = Index_Output_Variables[i];
    // Up to i=(SUB_OUTPUT_VARIABLES-1) are set through command line values:
    // -n [SUB_OUTPUT_VARIABLES] 
  
  /* Some implementations of this code require to alloc memmory according to
     a number of state variables that can change dynamically. In that case, this
     3 lines of code should be commented out. Notice alse the corresponding lines
     in the function below:
     void P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( Parameter_Table * T )
     Only if the number of dynamic state model variables never cange during execution,
     these lines of code make sense here
  */
  /* Initial Conditions: MODEL STATE VARIABLES */
  // for (i=0; i < MODEL_STATE_VARIABLES; i++){
  //    AssignLabel_to_Model_Variables(i, Table->Model_Variable_Name[i], Table);
  // }

  if(Table->TYPE_of_NETWORK == 0) {
    /// Setting up Constant Metapopulation Connectivity Matrix:
    for(a=0; a<Table->No_of_SPECIES; a++) 
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->No_of_CELLS; j++)
	  if (j != i) 
	    Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Mu;
	  else
	    Table->Metapop_Connectivity_Matrix[a][i][j] = 0.0; 
  }
  else {
    for(a=0; a<Table->No_of_SPECIES; a++) 
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->No_of_NEIGHBORS; j++)
	    Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Mu;
	  
  }
  /* END -------------------------------------------------*/
}
      

void P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( Parameter_Table * Table )
{
  int i, a;

  for (i=0; i < Table->MODEL_INPUT_PARAMETERS; i++){
    free( Table->Code_Parameters[i] );
    free( Table->Name_Parameters[i] );
    free( Table->Symbol_Parameters[i] ); 
    free( Table->Symbol_Parameters_Greek_CPGPLOT[i] ); 
  }
  free(Table->Code_Parameters);
  free(Table->Name_Parameters);
  free(Table->Symbol_Parameters); 
  free(Table->Symbol_Parameters_Greek_CPGPLOT); 
  
  free(Table->Default_Vector_Parameters);

  free(Table->Index);
  free(Table->Vector_Parameters);

  for (i=0; i < OUTPUT_VARIABLES_MAXIMUM; i++){
    free(Table->Output_Variable_Name[i]);
    free(Table->Output_Variable_Symbol[i]);
  }
  free(Table->Output_Variable_Symbol);
  free(Table->Output_Variable_Name);

  for (i=0; i<Table->MODEL_STATE_VARIABLES; i++){
     free(Table->Model_Variable_Name[i]);
     free(Table->Model_Variable_Symbol[i]);
  }
  free(Table->Model_Variable_Name);
  free(Table->Model_Variable_Symbol);

  /* Some implementations of this code require to alloc memmory according to
     a number of model variables that can change dynamically.

     For most calculations,  I opted for allocating/de-allocating this part of
     the data structure only if needed
     (see, for instance, ./MODEL_CALCULATIONS/TIME_EVO_DETERMINISTIC/MODEL.c)
     regardless the number of MODEL STATE VARIABLES changing or not during execution.
  */

  // for (i=0; i<MODEL_STATE_VARIABLES; i++){
  //   free(Table->Model_Variable_Name[i]);
  // }
  // free(Table->Model_Variable_Name);
  // free( Table->Vector_Model_Variables );
  // free( Table->Vector_Model_Variables_Stationarity );
  // free( Table->Vector_Model_Variables_MultiStability[0] );
  // free( Table->Vector_Model_Variables_MultiStability[1] );
  // free( Table->Vector_Model_Variables_MultiStability[2] );
  free(Table->Vector_Model_Variables_MultiStability);

  free(Table->OUTPUT_VARIABLE_INDEX);
  free(Table->Vector_Output_Variables);
  for (i=0; i<Table->SUB_OUTPUT_VARIABLES; i++) free(Table->Matrix_Output_Variables[i]);
  free(Table->Matrix_Output_Variables);

  free(Table->Default_Vector_Output_Variables);
 
  for(a=0; a<Table->No_of_SPECIES; a++) {  
    for(i=0; i<Table->No_of_CELLS; i++) 
      free(Table->Metapop_Connectivity_Matrix[a][i]); 
    
    free(Table->Metapop_Connectivity_Matrix[a]); 
  }
  free(Table->Metapop_Connectivity_Matrix); 
}

void Parameter_Table_Index_Update(int * Index, int N, Parameter_Table * P)
{
  int i;
  for(i=0; i<N; i++) P->Index[i] = Index[i];
}

/*
   The purpose of this simple function is just to upload
   ALL input parameters, controling both model definition
   and running simulations, which are defined as global variables,
   into the corresponding Parameter_Table Structure
*/
void Parameter_Values_into_Parameter_Table(Parameter_Table * P)
{
  int i_POP; 
  
  /* Setting up local populations */
  // for(i_POP = 0; i_POP<No_of_CELLS; i_POP++) {
  //  P->N_0[i_POP] = INITIAL_TOTAL_POPULATION;
  //  P->N_1[i_POP] = INITIAL_TOTAL_POPULATION;
  //  P->N_2[i_POP] = INITIAL_TOTAL_POPULATION;
  //  P->N_3[i_POP] = INITIAL_TOTAL_POPULATION;
  // }
  
  P->Mu                = Mu; 
  
  /* P->Beta     = Beta;     /\*  1 *\/     */
  /* P->Delta_0  = Delta_0;  /\*  2 *\/ */
  /* P->Alpha    = Alpha;    /\*  3 *\/      */
  /* P->Delta_1  = Delta_1;  /\*  4 *\/     */
  /* P->Gamma_1  = Gamma_1;  /\*  5 *\/ */
  /* P->Tau_1    = Tau_1;    /\*  6 *\/     */
  /* P->Y_P      = Y_P;      /\*  7 *\/ */
  /* P->Imm      = Imm;      /\*  8 *\/ */
  
  P->No_of_IC = No_of_IC;
  P->TYPE_of_INITIAL_CONDITION = TYPE_of_INITIAL_CONDITION;
  P->INITIAL_TOTAL_POPULATION  = INITIAL_TOTAL_POPULATION;

  P->RESCALING_INITIAL_TOTAL_POPULATION = RESCALING_INITIAL_TOTAL_POPULATION;

  P->TYPE_of_ERROR_FUNCTION = TYPE_of_ERROR_FUNCTION;
  P->No_of_ERROR_PARAMETERS = No_of_ERROR_PARAMETERS;

  P->Error_Parameter_0 = Error_Parameter_0;
  P->Error_Parameter_1 = Error_Parameter_1;

  P->Err_0 = Err_0;
  P->Err_1 = Err_1;
  P->Err_2 = Err_2;
  P->Err_3 = Err_3;
  P->Err_4 = Err_4;
  P->Err_5 = Err_5;
  P->Err_6 = Err_6;
  P->Err_7 = Err_7;
  P->Err_8 = Err_8;
  P->Err_9 = Err_9;
  P->Err_10 = Err_11;
  P->Err_11 = Err_11;
  P->Err_12 = Err_12;
  P->Err_13 = Err_13;
  P->Err_14 = Err_14;
  P->Err_15 = Err_15;
  P->Err_16 = Err_16;

  /* Definition of type of network */
  P->TYPE_of_NETWORK  = TYPE_of_NETWORK;
                               /* 0: Fully Connected Network
				  1: Square Grid with Von Neuman neighborhood
				  More network structures under construction 
			       */
  P->No_of_NEIGHBORS    = No_of_NEIGHBORS;

  P->No_of_SPECIES      = No_of_SPECIES; 
}
