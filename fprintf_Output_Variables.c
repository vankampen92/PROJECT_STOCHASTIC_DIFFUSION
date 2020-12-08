#include "./Include/MODEL.h"

void fprintf_Output_Variables(FILE *fp, Parameter_Table * P)
{
  int i,k;
  
  fprintf(fp, "\n");

  fprintf(fp, "Represented Output Variables:\n");
  fprintf(fp, "-----------------------------\n");
  
    for(i=0; i<P->SUB_OUTPUT_VARIABLES; i++){
      k = P->OUTPUT_VARIABLE_INDEX[i];
      fprintf(fp, "%d: Variable Name: %s\n", k, P->Output_Variable_Name[k]);
    }

}

