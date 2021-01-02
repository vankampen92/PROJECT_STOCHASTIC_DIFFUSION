#include "./Include/MODEL.h"

void AssignCPGPLOT_Symbol_to_Model_Parameters(int j, char * Label, Parameter_Table *P)
{
  /* Short Labels for Model Parameters */

  char * p;
  Label[0] = '\0';

  switch(j)
    {
    case  0: p=strcat(Label,"\\gm");    
      break;
    case  1: p=strcat(Label,"N");    
      break;
    case  2: p=strcat(Label,"M");       
      break; 
    case  3: p=strcat(Label,"M\\dx\\u");    
      break;
    case  4: p=strcat(Label,"M\\dy\\u");       
      break; 
    case  5: p=strcat(Label,"S\\dR\\u");       
      break; 
      
    case  6: p=strcat(Label, "\\gl\\u(R)\\d\\d0\\u");    
      break;
    case  7: p=strcat(Label, "\\gd\\u(R)\\d\\d0\\u");    
      break; 
    case  8: p=strcat(Label, "\\gl\\u(R)\\d\\d1\\u");        
      break;
    case  9: p=strcat(Label, "\\gd\\u(R)\\d\\d1\\u");    
      break; 
    case 10: p=strcat(Label, "K\\dR\\u");       
      break; 
    
    default:
      printf(".... INVALID PARAMETER KEY [key=%d]\n", j);
      printf(".... The permited correspondences are (0 to 12):\n");
      printf("\n");
      fprintf_Model_Parameters(stdout, P);
      exit(0);
    }
}


void AssignSymbol_to_Model_Parameters(int j, char * Label, Parameter_Table *P)
{
  /* Short Labels for Model Parameters */

  char * p;
  Label[0] = '\0';

  switch(j)
    {
    case  0: p=strcat(Label,"Mu");    
      break;
    case  1: p=strcat(Label,"N");    
      break;
    case  2: p=strcat(Label,"M");       
      break; 
    case  3: p=strcat(Label,"M_x");    
      break;
    case  4: p=strcat(Label,"M_y");       
      break; 
    case  5: p=strcat(Label,"S_R");       
      break; 

    case  6: p=strcat(Label,"Lambda_R_0");    
      break;
    case  7: p=strcat(Label,"Delta_R_0");       
      break; 
    case  8: p=strcat(Label,"Lambda_R_1");    
      break;
    case  9: p=strcat(Label,"Delta_R_1");       
      break; 
    case 10: p=strcat(Label,"K_R");       
      break; 
    
    default:
      printf(".... INVALID PARAMETER KEY [key=%d]\n", j);
      printf(".... The permited correspondences are (0 to 12)\n");
      printf("\n");
      fprintf_Model_Parameters(stdout, P);
      exit(0);
    }
}
