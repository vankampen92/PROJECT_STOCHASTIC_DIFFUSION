#include "./Include/MODEL.h"

extern int No_of_PARAMETERS; 

void fprintf_Model_Parameters(FILE * fp, Parameter_Table * P)
{
  int i,k;

  fprintf(fp, "Input Parameter:\t\t[-Input Code\tDefault Value]\n");
  fprintf(fp, "------------------------------------------------\n");

  for(i=0; i < No_of_PARAMETERS; i++){
      k = P->Index[i];
      fprintf(fp, "%d: %s:\t[ %s\t%g]\n", 
	      k, P->Name_Parameters[k], P->Code_Parameters[k], P->Default_Vector_Parameters[k] );
  }
  
}

