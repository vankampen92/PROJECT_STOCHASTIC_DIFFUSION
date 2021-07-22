#include <MODEL.h>

void AssignSymbol_to_Model_Variables(int j, char * Label, Parameter_Table * Table)
{
  char *pFile;
  int n,m; 
  char * pCell= (char *)calloc(10, sizeof(char) );

  n = j%Table->LOCAL_STATE_VARIABLES; 
  m = j/Table->LOCAL_STATE_VARIABLES;

  char * L[] = {"V", "R", "G"}; 
  
  Label[0]='\0';
  if( m < Table->No_of_CELLS ) { 
    
    sprintf(pCell, "%d", m);
    pFile = strcat(Label, "n\\d");
    pFile = strcat(Label, pCell);
    pFile = strcat(Label, "\\u");
    pFile = strcat(Label, "[");
   
    pFile = strcat(Label, L[n]);

    pFile = strcat(Label, "]");
    
  }  
  else{
    printf(".... INVALID VARIABLE KEY [key = %d]\n", j);
    printf(".... Provided Model Variable Codes have been correcly defined,\n");
    printf(".... the permited correspondences are:\n");
    printf(".... from key = 0 to key = %d\n", Table->MODEL_STATE_VARIABLES-1);
    printf(".... The program will exit\n");
    exit(0);
  }
  
  free(pCell); 
}

