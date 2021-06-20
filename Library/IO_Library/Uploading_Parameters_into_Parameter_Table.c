#include <MODEL.h>
#include "IO_Procedures_Aux.h"

void Uploading_Model_Parameters_into_Parameter_Table(Parameter_Table * Table, double ** Data,
						     int k)
{
  int i, key;
  int No_of_PARAMETERS; 
  
  No_of_PARAMETERS = Table->S->No_of_PARAMETERS;
  
  for(i=0; i < No_of_PARAMETERS; i++) {
    key = Table->S->Parameter_Index[i];
    AssignVectorEntry_to_Structure(Table , key, Data[k][i]);
  }
      
}

void Uploading_Full_Model_Parameters_into_Parameter_Table(Parameter_Table * Table, double ** Data,
							  int k)
{
  int n, i, key;
  int No_of_PARAMETERS; 
  
  No_of_PARAMETERS = Table->TOTAL_No_of_MODEL_PARAMETERS;  
  
  for(i=0; i < No_of_PARAMETERS; i++) {
    key = Table->Index[i];
    AssignVectorEntry_to_Structure(Table , key, Data[k][i]);
  }
  n = i; 
  
  if(Table->No_of_ERROR_PARAMETERS > 0) {
    for(i=0; i < Table->No_of_ERROR_PARAMETERS; i++) {
     key = Table->E_Space->Parameter_Index[i];
     Vector_Entry_into_Error_Model_Table(Data[k][n++], key, Table);
    }
  }

  if(Table->No_of_IC > 0) {
    for(i=0; i < Table->No_of_IC; i++) {
     key = Table->IC_Space->Parameter_Index[i];
     Vector_Entry_into_Model_Variable_Initial_Condition_Table(Data[k][n++], key, Table);
    }
  }
}


