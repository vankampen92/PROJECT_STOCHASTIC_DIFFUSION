#include <MODEL.h>

void AssignLabel_to_Output_Variables(int j, char * Label, Parameter_Table * P)
{
  char * p;
  Label[0] = '\0';

  int PUTA = OUTPUT_VARIABLES_GENUINE;
  printf("2: (j-%d) The number of output variables genuine is %d\n", j, PUTA);
 
  
  if (j >= OUTPUT_VARIABLES_GENUINE) {
    j -= OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    AssignLabel_to_Model_Variables(j, Label, P);
  }
  else {
    
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
      case  3:
	p = strcat(Label , "n(Sp_0)");   /*  0: Total No of Individals (Sp 0) */
        break;
      case  4:
	p = strcat(Label , "n(Sp_1)");   /*  1: Total No of Individals (Sp 1) */
        break;
      
      default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 4\n");
        exit(0);
      }
  }
}

void AssignLongLabel_to_Output_Variables(int j, char * Label, Parameter_Table * P)
{
  char * p;
  Label[0] = '\0';

  int PUTA = OUTPUT_VARIABLES_GENUINE;
  printf("3: (j-%d) The number of output variables genuine is %d\n", j, PUTA);
  
  if (j >= OUTPUT_VARIABLES_GENUINE) {
    j -= OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    AssignLabel_to_Model_Variables(j, Label, P);
  }
  else {
    
    switch(j)
      {
      case  0:
	p = strcat(Label , "Average Density");                   /*  0: S */
        break;
      case  1:
	p = strcat(Label , "Standard Deviation of the Density"); /*  0: S */
        break;
      case  2:
	p = strcat(Label , "Total No of Individuals");           /*  0: S */
        break;
      case  3:
	p = strcat(Label , "Total No of Individals (Sp 0)"); /*  0: S */
        break;
      case  4:
	p = strcat(Label , "Total No of Individals (Sp 1)");           /*  0: S */
        break;

      default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 4\n");
        exit(0);
      }
  }
}


