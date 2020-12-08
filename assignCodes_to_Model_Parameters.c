#include "./Include/MODEL.h"

void AssignCodes_to_Model_Parameters(int j, char * Label, Parameter_Table *P)
{
  char * p;
  Label[0] = '\0';

  switch(j)
    {
    case  0:  p = strcat(Label, "-Hu");  
	break;
    case  1:  p = strcat(Label, "-HN");  
      break;
    case  2:  p = strcat(Label, "-HM");
      break;
    case  3:  p = strcat(Label, "-HX");  
      break;
    case  4:  p = strcat(Label, "-HY");
      break;
    case  5:  p = strcat(Label, "-HS");
      break;
      
    default:
      printf(".... INVALID PARAMETER KEY [key = %d]\n", j);
      printf(".... The permited correspondences are 0 to 6:\n");
      printf("\n");
      fprintf_Model_Parameters(stdout, P);
      exit(0);
    }
}
