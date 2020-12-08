#include <MODEL.h>

void AssignLabel_to_Output_Variables(int j, char * Label, Parameter_Table * P)
{
  char * p;
  Label[0] = '\0';
  if (j >= OUTPUT_VARIABLES_GENUINE) {
    j -= OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    AssignLabel_to_Model_Variables(j, Label, P);
  }
  else {       /* SEI1I2AYR fraction and more  */
    switch(j)
      {
      case  0:  p = strcat(Label , "<n>");         /*  0: Density: No of Individuals per Cell */
        break;
      case  1:  p = strcat(Label , "s");         /*  1: SDV: No of Individuals per Cell */
        break;
      case  2:  p = strcat(Label , "N");         /*  2: Total No of Individuals */
        break;
      
      default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 2\n");
        exit(0);
      }
  }
}

void AssignLongLabel_to_Output_Variables(int j, char * Label, Parameter_Table * P)
{
  char * p;
  Label[0] = '\0';
  if (j >= OUTPUT_VARIABLES_GENUINE) {
    j -= OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    AssignLabel_to_Model_Variables(j, Label, P);
  }
  else {       /* SEI1I2AYR fraction and more  */
    switch(j)
      {
      case  0:  p = strcat(Label , "Average Density");                   /*  0: S */
        break;
      case  1:  p = strcat(Label , "Standard Deviation of the Density"); /*  0: S */
        break;
      case  2:  p = strcat(Label , "Total No of Individuals");           /*  0: S */
        break;

      default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 2\n");
        exit(0);
      }
  }
}


