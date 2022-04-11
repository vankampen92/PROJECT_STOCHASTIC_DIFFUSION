#include <MODEL.h>

double Model_Variable_Initial_Condition_into_Vector_Entry_Table(int j, Parameter_Table * Table)
{
  double Value;
  
  if( j < Table->MODEL_STATE_VARIABLES)  {
    Value = Table->Vector_Model_Variables_Time_0[j];
  }
  else{
    printf(".... INVALID VARIABLE KEY [key = %d]\n", j);
    printf(".... Provided Model Variable Codes have been correcly defined,\n");
    printf(".... the permited correspondences are:\n");
    printf(".... from key = 0 to key = %d\n", Table->MODEL_STATE_VARIABLES-1);
    printf(".... The program will exit\n");
    exit(0);
  }

  return(Value);
}

double Model_Variable_Initial_Condition_into_Vector_Entry_Model(int j, Parameter_Model * Table)
{
  double Value;

  if( j < Table->MODEL_STATE_VARIABLES)  {
    Value              = Table->Vector_Model_Variables_Time_0[j];
  }
  else{
    printf(".... INVALID VARIABLE KEY [key = %d]\n", j);
    printf(".... Provided Model Variable Codes have been correcly defined,\n");
    printf(".... the permited correspondences are:\n");
    printf(".... from key = 0 to key = %d\n", Table->MODEL_STATE_VARIABLES-1);
    printf(".... The program will exit\n");
    exit(0);
  }

  return(Value);
}
