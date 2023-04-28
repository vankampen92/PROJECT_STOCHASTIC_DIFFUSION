#include <MODEL.h>

void AssignLabel_to_Model_Variables(int j, char * Label, Parameter_Table * Table)
{
  char * pFile;
  int i,n,m; 
  char * pCell= (char *)calloc(10, sizeof(char) );
  char * L = (char *)calloc(10, sizeof(char) );
  
  n = j%Table->LOCAL_STATE_VARIABLES; /* Resource Type n  */ 
  m = j/Table->LOCAL_STATE_VARIABLES; /* Number of Cell m */

  L[0]='\0';
  if( n < Table->No_of_RESOURCES ) {
    sprintf(pCell, "%d", n);
    pFile = strcat(L, "AR_");
    pFile = strcat(L, pCell);
  }

  Label[0]='\0';
  if(m < Table->No_of_CELLS ) { 
    
    sprintf(pCell, "%d", m);
    pFile = strcat(Label, "n_");
    pFile = strcat(Label, pCell);

    pFile = strcat(Label, "[");
    pFile = strcat(Label, L);
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

  free(pCell); free(L);
}
