#include <MODEL.h>

void AssignSymbol_to_Output_Variables(int j, char * Label, Parameter_Table * Table)
{
  char * p;
  Label[0] = '\0';
  if (j >= OUTPUT_VARIABLES_GENUINE) {
    j -= OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    AssignSymbol_to_Model_Variables(j, Label, Table);
  }
  else {        /* SEI1I2AYR fraction and more  */
    switch(j)
      {
      case  0:  p = strcat(Label , "x");        /*  0 */
        break;
      case  1:  p = strcat(Label , "s");         /*  1: SDV: No of Individuals per Cell */
        break;
      case  2:  p = strcat(Label , "N");         /*  2: Total No of Individuals */
        break;
	
      default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 13\n");
        exit(0);
      }
  }
}

void AssignCPGPLOT_Symbol_to_Output_Variables(int j, char * Label, Parameter_Table * Table)
{
  char * p;
  Label[0] = '\0';
  if (j >= OUTPUT_VARIABLES_GENUINE) {
    j -= OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    AssignSymbol_to_Model_Variables(j, Label, Table);
  }
  else {
    switch(j)
      {
      case  0:  p = strcat(Label , "x");       /*  0 */
	break;
      case  1:  p = strcat(Label , "s");         /*  1: SDV: No of Individuals per Cell */
        break;
      case  2:  p = strcat(Label , "N");         /*  2: Total No of Individuals */
        break;
	
      default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 13\n");
        exit(0);
      }
  }
}
