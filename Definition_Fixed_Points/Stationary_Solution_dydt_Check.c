#include <MODEL.h>

void Stationary_Solution_dvdt_Check (Parameter_Table * Table )
{
  /* Only valid on one single patch */
  int i, k;
  int Positive; 
  char * Label; 
  double value;
  
  assert(Table->No_of_CELLS == 1 && Table->MODEL_STATE_VARIABLES == 4); 
  
  Label = (char *)calloc( 20, sizeof(char) ); 
    
  Positive = 1; 
  for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) 
    if (Table->Vector_Model_Variables_Stationarity[i] < 0.0) Positive = 0; 

  if( Positive == 0 ) {
    printf("Model Variables at stationarity: feasible (postive) stationarity point\n");
    
    double * dydt = (double *)calloc(Table->MODEL_STATE_VARIABLES, sizeof(double));
    
    void * params = Table; 
    int SIGNAL    = function (0.0, Table->Vector_Model_Variables_Stationarity, dydt, void *params);

    assert(SIGNAL == GSL_SUCCESS);

    for (k=0; k < Table->MODEL_STATE_VARIABLES; k++) {
      AssignLabel_to_Model_Variables(k, Label, Table);
      printf("y[%s]    = %g\t", Label, Table->Vector_Model_Variables_Stationarity[k]);
      printf("dydt[%s] = %g\n", Label, dydt[k]);
    }
    printf("\n");

    free(dydt);
  }
  else {
    
    printf("Model Variables at stationarity are not feasible!!!\n");
  }

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