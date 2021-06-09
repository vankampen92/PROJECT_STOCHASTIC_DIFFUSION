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

