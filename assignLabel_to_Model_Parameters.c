#include "./Include/MODEL.h"

void AssignLabel_to_Model_Parameters(int j, char * Label, Parameter_Table *P)
{
  char * p;
  Label[0] = '\0';

  switch(j)
    {
    case  0:  p = strcat(Label, "Jumping Rate");  
	break;
    case  1:  p = strcat(Label, "No of INDIVIDUALS");  
      break;
    case  2:  p = strcat(Label, "No of CELLS");  
      break;
    case  3:  p = strcat(Label, "No of CELLS (in horizontal dimension)");  
      break;
    case  4:  p = strcat(Label, "No of CELLS (in vertical dimension)");  
      break;
    case  5:  p = strcat(Label, "No of (Resource) SPECIES");  /* Number of Resource Species */
      break;
    case  6:  p = strcat(Label, "External Immigration Rate (0)");  
      break;
    case  7:  p = strcat(Label, "Decaying Rate (0)");
      break;
    case  8:  p = strcat(Label, "External Immigration Rate (1)");  
      break;
    case  9:  p = strcat(Label, "Decaying Rate (1)");
      break;
    case 10:  p = strcat(Label, "Patch Carrying Capacity (Resources)");   /* Patch Carrying Capacity */
      break;
    
      
    default:
      printf(".... INVALID PARAMETER KEY [key = %d]\n", j);
      printf(".... The permited correspondences are (0 to 6):\n");
      printf("\n");
      fprintf_Model_Parameters(stdout, P);
      exit(0);
    }
}
