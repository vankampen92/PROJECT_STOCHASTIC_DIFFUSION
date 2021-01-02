#include <MODEL.h>

void AssignLabel_to_Output_Variables(int j, char * Label, Parameter_Table * Table)
{
  char * p;
  Label[0] = '\0';
  
  if (j >= Table->OUTPUT_VARIABLES_GENUINE) {
    j -= Table->OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    AssignLabel_to_Model_Variables(j, Label, Table);
  }
  else if (j < Table->No_of_RESOURCES ) {
    char * n_Label = (char *)calloc(10, sizeof(char) );
    p = strcat(Label, "n_R[");
    n_Label[0]='\0';
    sprintf(n_Label, "%d", j);
    p = strcat(Label, n_Label);
    p = strcat(Label, "]");
    free(n_Label);   
  }
  else if (j < Table->OUTPUT_VARIABLES_GENUINE) {
    j -= Table->No_of_RESOURCES;
    
    switch(j)
      {
      case  0:
	p = strcat(Label , "<n>");         /*  0: Density: No of Individuals per Cell */
        break;
      case  1:
	p = strcat(Label , "s");         /*  1: SDV: No of Individuals per Cell */
        break;
      case  2:
	p = strcat(Label , "N");         /*  2: Total No of Individuals */
        break;
      
      default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 4\n");
        exit(0);
      }
  }
}

void AssignLongLabel_to_Output_Variables(int j, char * Label, Parameter_Table * Table)
{
  char * p;
  Label[0] = '\0';


  if (j >= Table->OUTPUT_VARIABLES_GENUINE) {
    j -= Table->OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    AssignLabel_to_Model_Variables(j, Label, Table);
  }
  else if (j < Table->No_of_RESOURCES ) {
    char * n_Label = (char *)calloc(10, sizeof(char) );
    p = strcat(Label, "n_R[");
    n_Label[0]='\0';
    sprintf(n_Label, "%d", j);
    p = strcat(Label, n_Label);
    p = strcat(Label, "]");
    free(n_Label);   
  }
  else if (j < Table->OUTPUT_VARIABLES_GENUINE) {
    j -= Table->No_of_RESOURCES;
    
    switch(j) {

    case  0:
      p = strcat(Label , "Average Density");                   /*  0: S */
      break;
    case  1:
      p = strcat(Label , "Standard Deviation of the Density"); /*  0: S */
      break;
    case  2:
      p = strcat(Label , "Total No of Individuals");           /*  0: S */
      break;
      
    default:
      printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
      printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 4\n");
        exit(0);
    }
  }
}


