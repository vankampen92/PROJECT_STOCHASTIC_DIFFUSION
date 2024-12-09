#include <MODEL.h>

void AssignLabel_to_Model_Variables(int j, char * Label, Parameter_Table * Table)
{
  char * pFile;
  int i,k, n,m; 
  char * Lo_Va = (char *)calloc(10, sizeof(char) );
  char * pCell = (char *)calloc(10, sizeof(char) );
  
  n = j%Table->LOCAL_STATE_VARIABLES; 
  m = j/Table->LOCAL_STATE_VARIABLES;  

  assert(Table->LOCAL_STATE_VARIABLES == Table->No_of_RESOURCES);

  Calcualte_Strain_and_Profile(Table, n, &i, &k);

  Lo_Va[0]='\0';
  sprintf(Lo_Va, "%d ,", i);  /* ith Strain Type */
  pFile = strcat(Lo_Va, "[");
  sprintf(Lo_Va, "%d", k);    /* kth Profile Type (within the i-th Strain Type */
  pFile = strcat(Lo_Va, "]");
  
  Label[0]='\0';
  if( m < Table->No_of_CELLS ) { 
    
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

