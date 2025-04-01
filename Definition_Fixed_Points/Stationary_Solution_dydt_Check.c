#include <MODEL.h>

void Stationary_Solution_dvdt_Check (Parameter_Table * Table )
{
  /* Only valid on one single patch */
  int i, k;
  int Positive; 
  char * Label; 
  double value;
  void * params = (void *)Table;
    
  /* Ensure the model is only dealing with one cell */
  assert(Table->No_of_CELLS == 1); 

  /* Allocate memory for dydt */
  double * dydt = (double *)calloc(Table->MODEL_STATE_VARIABLES, sizeof(double));
  
  /* Allocate memory for the label */
  Label = (char *)calloc( 20, sizeof(char) ); 
    
  Positive = 1; 
  /* Check if all model variables at stationarity are positive */
  for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) 
    if (Table->Vector_Model_Variables_Stationarity[i] < 0.0) Positive = 0; 

  if( Positive == 1 ) {
    printf("Model Variables at stationarity: feasible (positive) stationarity point\n");

    /* Evaluate the function at stationarity */
    int SIGNAL    = function (0.0, Table->Vector_Model_Variables_Stationarity, dydt, params);

    /* Ensure the function evaluation was successful */
    assert(SIGNAL == GSL_SUCCESS);

    /* Print the model variables and their derivatives at stationarity */
    for (k=0; k < Table->MODEL_STATE_VARIABLES; k++) {
      AssignLabel_to_Model_Variables(k, Label, Table);
      printf("y[%s]    = %g\t", Label, Table->Vector_Model_Variables_Stationarity[k]);
      printf("dydt[%s] = %g\n", Label, dydt[k]);
    }
    printf("\n");
  }
  else {
    /* Print the model variables at stationarity if they are not feasible */
    for (k=0; k < Table->MODEL_STATE_VARIABLES; k++) {
      AssignLabel_to_Model_Variables(k, Label, Table);
      printf("y[%s]    = %g\t", Label, Table->Vector_Model_Variables_Stationarity[k]);
    }
    printf("\n");    
    printf("Model Variables at stationarity are not feasible!!!\n");
  }
  /* Free the allocated memory for Label */
  free(Label);  
  
  printf("Model Parameters:\n");
  /* Allocate memory for the label again */
  Label = (char *)calloc( 20, sizeof(char) );
  /* Print all model parameters */
  for (k=0; k < Table->TOTAL_No_of_MODEL_PARAMETERS; k++) {
      i = Table->Index[k]; 
      AssignSymbol_to_Model_Parameters(i, Label, Table);
      value = AssignStructValue_to_VectorEntry(i, Table);
      printf("%s = %g\n", Label, value);
  }
  printf("\n");
  /* Free the allocated memory for Label */
  free(Label);
  /* Free the allocated memory for dydt */
  free(dydt);    
}