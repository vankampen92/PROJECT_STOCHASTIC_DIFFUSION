#include "./Include/MODEL.h"

void AssignLabel_to_Model_Parameters__LATEX(int j, char * Label, Parameter_Table *P)
{
  char * p;
  Label[0] = '\0';

  switch(j)
    {
    
    case  0:  p = strcat(Label , "Diffusion Jumping Rate");    /*  0 */
      break; 
    case  1:  p = strcat(Label , "No of INDIVIDUALS");         /*  1 */
      break;
    case  2:  p = strcat(Label , "No_of_CELLS");               /*  2 */
      break;
    case  3:  p = strcat(Label , "No of CELSS (X Dimension)"); /*  3 */
      break;
    case  4:  p = strcat(Label , "No_of_CELLS (Y Dimension)"); /*  4 */
      break;
    case  5:  p = strcat(Label , "No_of_RESOURCES"); /*  4 */
      break;
    case  6:  p = strcat(Label, "External Immigration Rate (0)");  
      break;
    case  7:  p = strcat(Label, "Decaying Rate (0)");
      break;
    case  8:  p = strcat(Label, "External Immigration Rate (1)");  
      break;
    case  9:  p = strcat(Label, "Decaying Rate (1)");
      break;
    case 10:  p = strcat(Label, "Patch Carrying Capacity");   /* Patch Carrying Capacity */
      break;
    case 11: p=strcat(Label,  "Resource Local Reproduction Rate");        /* -H5 */
      break;

    case 12: p=strcat(Label,  "Consumer External Immigration Rate (0)");     /* -H5 */
      break;
    case 13: p=strcat(Label,  "Consumer Death Rate (0)");      /* -H6 */
      break;
    case 14: p=strcat(Label,  "Consumer External Immigration Rate (1)");     /* -H7 */
      break;
    case 15: p=strcat(Label,  "Consumer Death Rate (1)");      /* -H8 */ 
      break; 
      
    case 16: p=strcat(Label,  "Comsumer Attack Rate (0)");      /* -H9 */
      break;
    case 17: p=strcat(Label,  "Nu = 1/Tau (One over the handling time)");  /* -H10 */
      break;

    case 18: p=strcat(Label,  "Triplet formation rate (0)");        /* -H11 */
      break;
    case 19: p=strcat(Label,  "Triplet desintegration (0)");        /* -H12 */
      break;

    case 20: p=strcat(Label,  "Consumer Movement Rate (0)");        /* -H13 */
      break;
    
      
    default:
      printf(".... INVALID PARAMETER KEY [key = %d]\n", j);
      printf(".... The permited correspondences are:\n");
      printf("\n");
      fprintf_Model_Parameters(stdout, P);
      exit(0);
    }
}

void AssignLabel_to_Model_Parameters__LATEX__SYMBOL(int j, char * Label, Parameter_Table *P)
{
  char * p;
  Label[0] = '\0';

  switch(j)
    {

    case  0:  p = strcat(Label, "$\\mu$");
      break;
    case  1:  p = strcat(Label, "$N$");      
      break;
    case  2:  p = strcat(Label, "$M$");        
      break;
    case  3:  p = strcat(Label, "$M_X$");      
      break;
    case  4:  p = strcat(Label, "$M_Y$");        
      break;
    case  5:  p = strcat(Label, "$S_R$");        
      break;
    case  6:  p = strcat(Label, "$\\lambda^{(R)}_0$");  
      break;
    case  7:  p = strcat(Label, "$\\delta^{(R)}_0$");
      break;
    case  8:  p = strcat(Label, "$\\lambda^{(R)}_1$");  
      break;
    case  9:  p = strcat(Label, "$\\delta^{(R)}_1$");
      break;
    case 10:  p = strcat(Label, "$K_R$");   /* Patch Carrying Capacity */
      break;
    case 11: p=strcat(Label,"$\\beta_R$");        /* -H5 */
      break;

    case 12: p=strcat(Label,"$\\lambda^{(C)}_0$");     /* -H5 */
      break;
    case 13: p=strcat(Label,"$\\delta^{(C)}_0$");     /* -H6 */
      break;
    case 14: p=strcat(Label,"$\\lambda^{(C)}_1$");    /* -H7 */
      break;
    case 15: p=strcat(Label,"$\\delta^{(C)}_1$");     /* -H8 */ 
      break; 
      
    case 16: p=strcat(Label,"$\\alpha^{(C)}_0$");     /* -H9 */
      break;
    case 17: p=strcat(Label,"$\\nu^{(C)}_0$");        /* -H10 */
      break;

    case 18: p=strcat(Label,"$\\chi^{(C)}_0$");        /* -H11 */
      break;
    case 19: p=strcat(Label,"$\\eta^{(C)}_0$");        /* -H12 */
      break; 

    case 20: p=strcat(Label,"$\\mu^{(C)}$");           /* -H13 */
      break; 
   
    default:
      printf(".... INVALID PARAMETER KEY [key = %d]\n", j);
      printf(".... The permited correspondences are:\n");
      printf("\n");
      fprintf_Model_Parameters(stdout, P);
      exit(0);
    }
}
