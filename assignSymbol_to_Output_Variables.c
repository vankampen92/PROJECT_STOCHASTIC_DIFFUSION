#include <MODEL.h>

void AssignSymbol_to_Output_Variables(int j, char * Label, Parameter_Table * Table)
{
  char * p;
  Label[0] = '\0';

  int PUTA = OUTPUT_VARIABLES_GENUINE;
  printf("0: (j-%d) The number of output variables genuine is %d\n", j, PUTA);
  
  if (j >= OUTPUT_VARIABLES_GENUINE) {
    j -= OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    AssignSymbol_to_Model_Variables(j, Label, Table);
  }
  else {     
    switch(j)
      {
      case  0:
	p = strcat(Label , "x");        /*  0 */
        break;
      case  1:
	p = strcat(Label , "s");         /*  1: SDV: No of Individuals per Cell */
        break;
      case  2:
	p = strcat(Label , "N");         /*  2: Total No of Individuals */
        break;
      case  3:
	p = strcat(Label , "n(Sp_0)");   /*  3: Total No of Individals (Sp 0) */
        break;
      case  4:
	p = strcat(Label , "n(Sp_1)");   /*  4: Total No of Individals (Sp 1) */
        break;
	
      default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 4\n");
        exit(0);
      }
  }
}

void AssignCPGPLOT_Symbol_to_Output_Variables(int j, char * Label, Parameter_Table * Table)
{
  char * p;
  Label[0] = '\0';

  int PUTA = OUTPUT_VARIABLES_GENUINE;
  printf("1: (j-%d) The number of output variables genuine is %d\n", j, PUTA);
  
  if (j >= OUTPUT_VARIABLES_GENUINE) {
    j -= OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    AssignSymbol_to_Model_Variables(j, Label, Table);
  }
  else {
    
    switch(j)
      {
      case  0:
	p = strcat(Label , "x");       /*  0 */
	break;
      case  1:
	p = strcat(Label , "s");         /*  1: SDV: No of Individuals per Cell */
        break;
      case  2:
	p = strcat(Label , "N");         /*  2: Total No of Individuals */
        break;
      case  3:
	p = strcat(Label , "N[0]");      /*  3: Total No of Individals (Sp 0) */
        break;
      case  4:
	p = strcat(Label , "N[1]");      /*  4: Total No of Individals (Sp 1) */
        break;
	
      default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 4\n");
        exit(0);
      }
  }
}
