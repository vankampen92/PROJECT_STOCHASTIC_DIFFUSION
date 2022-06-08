#include <MODEL.h>

void Stationary_Solution_Feasibility_Control (Parameter_Table * Table )
{
  /* Only valid on one single patch */
  int i, k;
  int Positive; 
  char * Label; 
  double value;
  
  assert(Table->No_of_CELLS == 1); 
  
  Positive = 1; 
  for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) 
    if (Table->Vector_Model_Variables_Stationarity[i] < 0.0) Positive = 0; 

  if( Positive == 0 ) {
    printf("Model Variables at stationarity\n");
    Label = (char *)calloc( 20, sizeof(char) );
    for (k=0; k < Table->MODEL_STATE_VARIABLES; k++) {
      AssignLabel_to_Model_Variables(k, Label, Table);
      printf("y[%s] = %g\n", Label, Table->Vector_Model_Variables_Stationarity[k]);
    }
    printf("\n");
    free(Label);
    
    printf("Model Parameters:\n");
    Label = (char *)calloc( 20, sizeof(char) );
    for (k=0; k < Table->TOTAL_No_of_MODEL_PARAMETERS; k++) {
      i = Table->Index[k]; 
      AssignSymbol_to_Model_Parameters(i, Label, Table);
      value = AssignStructValue_to_VectorEntry(i, Table);
      printf("%s = %g\n", Label, value);
    }
    printf("\n");
    free(Label);
  }
 
  assert(Positive == 1);
}

void assert_positive_model_parameters (Parameter_Table * Table )
{
  /* Only valid on one single patch */
  int i, k;
  int Positive; 
  char * Label; 
  double value;
  
  assert(Table->No_of_CELLS == 1); 
  
  Positive = 1;
  
  printf("Model Parameters:\n");
  Label = (char *)calloc( 20, sizeof(char) );
  for (k=0; k < Table->TOTAL_No_of_MODEL_PARAMETERS; k++) {
    i = Table->Index[k]; 
    AssignSymbol_to_Model_Parameters(i, Label, Table);
    value = AssignStructValue_to_VectorEntry(i, Table);
    printf("%s = %g\n", Label, value);
    if( value < 0.0) Positive = 0; 
  }
  printf("\n");
  free(Label);
  
  assert(Positive == 1);
}
