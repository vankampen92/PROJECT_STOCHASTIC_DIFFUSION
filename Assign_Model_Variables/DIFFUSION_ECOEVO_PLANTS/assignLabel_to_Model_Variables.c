#include <MODEL.h>

void AssignLabel_to_Model_Variables(int j, char * Label, Parameter_Table * Table)
{
  char * pFile;
  int i,k, n,m; 
  char * Lo_Va = (char *)calloc(10, sizeof(char) );
  char * pCell = (char *)calloc(10, sizeof(char) );
  
  n = j%Table->LOCAL_STATE_VARIABLES; 
  m = j/Table->LOCAL_STATE_VARIABLES;  

  if( n >= 2 * Table->No_of_RESOURCES ) {
    printf(".... INVALID VARIABLE KEY [key = %d]\n", j);
    printf(".... LOCAL VARIABLER KEY  [key = %d] (from 0 to %d-1)\n", n, 2*Table->No_of_RESOURCES);
    printf(".... Provided Model Variable Codes have been correcly defined,\n");
    printf(".... the permited correspondences are:\n");
    printf(".... from key = 0 to key = %d\n", Table->MODEL_STATE_VARIABLES-1);
    printf(".... The program will exit\n");
    exit(0);
  } 

  char * L[] = {"RP", "R"};

  i = n%2;
  k = n/2;

  Lo_Va[0]='\0';
  sprintf(Lo_Va, "%s", L[i]);  /* Stage: R or RP  */
  pFile = strcat(Lo_Va, "[");
  sprintf(Lo_Va, "%d", k);     /* Species or type */
  pFile = strcat(Lo_Va, "]");
  
  if(m < Table->No_of_CELLS ) { 
    
    pCell[0]='\n';
    sprintf(pCell, "%d", m);
    pFile = strcat(Label, "n_");
    pFile = strcat(Label, pCell);
    pFile = strcat(Label, "[");
    pFile = strcat(Label, Lo_Va);
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
  free(Lo_Va);
}
