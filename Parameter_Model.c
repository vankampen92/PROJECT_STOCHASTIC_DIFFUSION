#include "./Include/MODEL.h"

void Parameter_Model_Copy (Parameter_Model * P_Destination, Parameter_Model * P_Source)
{

  P_Destination->Mu         = P_Source->Mu;      /*  0 */ 
  P_Destination->Lambda_R_0 = P_Source->Lambda_R_0; 
  P_Destination->Delta_R_0  = P_Source->Delta_R_0; 
  P_Destination->Lambda_R_1 = P_Source->Lambda_R_1; 
  P_Destination->Delta_R_1 = P_Source->Delta_R_1;
  P_Destination->K_R        = P_Source->K_R; 
  
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
  
  P_Destination->MODEL_OUTPUT_VARIABLES = P_Source->MODEL_OUTPUT_VARIABLES;   /* Actual no of MODEL (output) VARIABLES */
  P_Destination->MODEL_INPUT_PARAMETERS = P_Source->MODEL_INPUT_PARAMETERS;  /* Actual no of MODEL (input) PARAMETERS */
  P_Destination->MODEL_STATE_VARIABLES  = P_Source->MODEL_STATE_VARIABLES;  /* Actual no of MODEL (state) VARIABLES  */

  P_Destination->TYPE_of_MODEL        = P_Source->TYPE_of_MODEL;
  P_Destination->Growth_Function_Type = P_Source->Growth_Function_Type;

  //P_Destination->Time_Scale_Unit      = Table->T->Time_Scale_Unit;

  P_Destination->No_of_CELLS       = P_Source->No_of_CELLS; 
  P_Destination->No_of_INDIVIDUALS = P_Source->No_of_INDIVIDUALS;
  P_Destination->No_of_CELLS_X     = P_Source->No_of_CELLS_X;
  P_Destination->No_of_CELLS_Y     = P_Source->No_of_CELLS_Y;
  
  /* Ex: 8 (without acummulating variables) */
  P_Destination->No_of_NEIGHBORS            = P_Source->No_of_NEIGHBORS;
  P_Destination->TYPE_of_NETWORK            = P_Source->TYPE_of_NETWORK;
  P_Destination->TOTAL_No_of_EVENTS         = P_Source->TOTAL_No_of_EVENTS;
  /* Total Number of Events within a patch, i.e., 25 * 4 = 100 */
  P_Destination->No_of_EVENTS               = P_Source->No_of_EVENTS;
  /* Number of Events within an age class, i.e., 25           */

  P_Destination->Metapop_Connectivity_Matrix = P_Source->Metapop_Connectivity_Matrix;

  P_Destination->No_of_RESOURCES = P_Source->No_of_RESOURCES; 
}

void  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N ( Parameter_Table * Table,
							Parameter_Model * P )
{
  /* This function transfer a subset of table parameters
     into the Parameter_Model structure. Parameter_Model parameters
     control the dynamical model.
  */

  P->Mu         = Table->Mu;       /*  0 */ 
  P->Lambda_R_0 = Table->Lambda_R_0; 
  P->Delta_R_0  = Table->Delta_R_0; 
  P->Lambda_R_1 = Table->Lambda_R_1; 
  P->Delta_R_1  = Table->Delta_R_1;
  P->K_R        = Table->K_R; 
  
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

  P->MODEL_OUTPUT_VARIABLES = Table->MODEL_OUTPUT_VARIABLES;   /* Actual no of MODEL (output) VARIABLES */
  P->MODEL_INPUT_PARAMETERS = Table->MODEL_INPUT_PARAMETERS;  /* Actual no of MODEL (input) PARAMETERS */
  P->MODEL_STATE_VARIABLES  = Table->MODEL_STATE_VARIABLES;  /* Actual no of MODEL (state) VARIABLES  */

  P->TYPE_of_MODEL        = Table->TYPE_of_MODEL;
  P->Growth_Function_Type = Table->Growth_Function_Type;

  //P->Time_Scale_Unit      = Table->T->Time_Scale_Unit;

  P->No_of_CELLS = Table->No_of_CELLS;
  P->No_of_CELLS       = Table->No_of_CELLS; 
  P->No_of_INDIVIDUALS = Table->No_of_INDIVIDUALS;
  P->No_of_CELLS_X     = Table->No_of_CELLS_X;
  P->No_of_CELLS_Y     = Table->No_of_CELLS_Y;
  
  P->No_of_NEIGHBORS            = Table->No_of_NEIGHBORS;
  P->TYPE_of_NETWORK            = Table->TYPE_of_NETWORK;
  P->TOTAL_No_of_EVENTS         = Table->TOTAL_No_of_EVENTS;
  /* Total Number of Events within a patch, i.e., 18 * 4 = 72 */
  P->No_of_EVENTS               = Table->No_of_EVENTS;
  /* Number of Events within an age class, i.e., 18           */

  P->Metapop_Connectivity_Matrix = Table->Metapop_Connectivity_Matrix;

  P->No_of_RESOURCES = Table->No_of_RESOURCES; 
}

void Parameter_Model_Copy_into_Parameter_Table (Parameter_Table * P_Destination, Parameter_Model * P_Source)
{
  P_Destination->Mu      = P_Source->Mu;       /*  0 */ 

  P_Destination->Lambda_R_0 = P_Source->Lambda_R_0; 
  P_Destination->Delta_R_0  = P_Source->Delta_R_0; 
  P_Destination->Lambda_R_1 = P_Source->Lambda_R_1; 
  P_Destination->Delta_R_1  = P_Source->Delta_R_1;
  P_Destination->K_R        = P_Source->K_R; 
  
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
  
  P_Destination->MODEL_OUTPUT_VARIABLES = P_Source->MODEL_OUTPUT_VARIABLES;   /* Actual no of MODEL (output) VARIABLES */
  P_Destination->MODEL_INPUT_PARAMETERS = P_Source->MODEL_INPUT_PARAMETERS;  /* Actual no of MODEL (input) PARAMETERS */
  P_Destination->MODEL_STATE_VARIABLES  = P_Source->MODEL_STATE_VARIABLES;  /* Actual no of MODEL (state) VARIABLES  */

  P_Destination->TYPE_of_MODEL        = P_Source->TYPE_of_MODEL;
  P_Destination->Growth_Function_Type = P_Source->Growth_Function_Type;

  //P->Time_Scale_Unit      = Table->T->Time_Scale_Unit;
  
  P_Destination->No_of_CELLS       = P_Source->No_of_CELLS; 
  P_Destination->No_of_INDIVIDUALS = P_Source->No_of_INDIVIDUALS;
  P_Destination->No_of_CELLS_X     = P_Source->No_of_CELLS_X;
  P_Destination->No_of_CELLS_Y     = P_Source->No_of_CELLS_Y;

  /* Ex: 8 (without acummulating variables) */
  P_Destination->No_of_NEIGHBORS            = P_Source->No_of_NEIGHBORS;
  P_Destination->TYPE_of_NETWORK            = P_Source->TYPE_of_NETWORK;
  P_Destination->TOTAL_No_of_EVENTS         = P_Source->TOTAL_No_of_EVENTS;
  /* Total Number of Events within a patch, i.e., 18 * 4 = 72 */
  P_Destination->No_of_EVENTS               = P_Source->No_of_EVENTS;
  /* Number of Events within an age class, i.e., 18           */

  P_Destination->Metapop_Connectivity_Matrix = P_Source->Metapop_Connectivity_Matrix;

  P_Destination->No_of_RESOURCES = P_Source->No_of_RESOURCES; 
}

void Parameter_Table_Copy_into_Parameter_Model (Parameter_Model * P_Destination, Parameter_Table * P_Source)
{

  P_Destination->Mu      = P_Source->Mu;       /*  0 */ 

  P_Destination->Lambda_R_0 = P_Source->Lambda_R_0; 
  P_Destination->Delta_R_0  = P_Source->Delta_R_0; 
  P_Destination->Lambda_R_1 = P_Source->Lambda_R_1; 
  P_Destination->Delta_R_1  = P_Source->Delta_R_1;
  P_Destination->K_R        = P_Source->K_R; 

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

  P_Destination->MODEL_OUTPUT_VARIABLES = P_Source->MODEL_OUTPUT_VARIABLES;   /* Actual no of MODEL (output) VARIABLES */
  P_Destination->MODEL_INPUT_PARAMETERS = P_Source->MODEL_INPUT_PARAMETERS;  /* Actual no of MODEL (input) PARAMETERS */
  P_Destination->MODEL_STATE_VARIABLES  = P_Source->MODEL_STATE_VARIABLES;  /* Actual no of MODEL (state) VARIABLES  */

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
  /* Total Number of Events within a patch, i.e., 18 * 4 = 72 */
  P_Destination->No_of_EVENTS               = P_Source->No_of_EVENTS;
  /* Number of Events within an age class, i.e., 18           */
  
  P_Destination->Metapop_Connectivity_Matrix = P_Source->Metapop_Connectivity_Matrix;

  P_Destination->No_of_RESOURCES = P_Source->No_of_RESOURCES; 
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
  case 10: value = (double)P->K_R;              /* Resource Carrying Capacity */ 
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
