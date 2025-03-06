#include "./Include/MODEL.h"

void AssignCodes_to_Model_Parameters(int j, char * Label, Parameter_Table *P)
{
  char * p;
  Label[0] = '\0';

  switch(j)
    {
    case  0:  p = strcat(Label, "-Hu");   /* Mu Per Capita Movement Rate */ /* Nu o be Nu_0 */  
	   break;
    case  1:  p = strcat(Label, "-HN");  
      break;
    case  2:  p = strcat(Label, "-HM");
      break;
    case  3:  p = strcat(Label, "-HX");  
      break;
    case  4:  p = strcat(Label, "-HY");
      break;
    case  5:  p = strcat(Label, "-HS");   /* Number of Resource Species */
      break;
    case  6:  p = strcat(Label, "-H0");   /* Lambda_R_0 */ 
      break;
    case  7:  p = strcat(Label, "-H1");   /* Delta_R_0 */
      break; 
    case  8:  p = strcat(Label, "-H2");   /* Lambda_R_1 */  
      break;
    case  9:  p = strcat(Label, "-H3");   /* Delta_R_1 */
      break;
    case 10:  p = strcat(Label, "-HK");   /* Patch Carrying Capacity */
      break;
    case 11:  p = strcat(Label, "-H4");   /* Beta_R */
      break;

    case 12:  p = strcat(Label, "-H5");   /* Lambda_C_0 */ 
      break;
    case 13:  p = strcat(Label, "-H6");    /* Delta_C_0 */   
      break;
#ifdef DIFFUSION_AZTECA_4D_1
    case 14:  p = strcat(Label, "-HuQ");   /* Lambda_C_1 */ /* Overloaded!!! */     
      break;
#else
    case 14:  p = strcat(Label, "-H7");    /* Lambda_C_1 */     
      break;
#endif
    case 15:  p = strcat(Label, "-H8");    /* Delta_C_1 */   
      break; 
      
    case 16:  p = strcat(Label, "-H9");    /* Alpha_C_0 */
      break;
    
    case 17:  p = strcat(Label, "-H10");   /* Nu_C_0    */
      break;

    case 18:  p = strcat(Label, "-H11");   /* Xhi_C_0 */
      break;
      
    case 19:  p = strcat(Label, "-H12");   /* Eta_C_0 */
      break;

    case 20:  p = strcat(Label, "-H13");   /* Mu_C */ /* Nu_1 */  
      break;

    case 21:  p = strcat(Label, "-H14");   
      break;
    case 22:  p = strcat(Label, "-H15");   
      break;
    case 23:  p = strcat(Label, "-H16");   
      break; 
    case 24:  p = strcat(Label, "-H17");   
      break;
    case 25:  p = strcat(Label, "-H18");   
      break;
    case 26:  p = strcat(Label, "-H19");   
      break;
    case 27:  p = strcat(Label, "-Hp1");   
      break; 
    case 28:  p = strcat(Label, "-Hp2");   
      break;
    case 29:  p = strcat(Label, "-H20");   
      break;

    default:
      printf(".... INVALID PARAMETER KEY [key = %d]\n", j);
      printf(".... The permited correspondences are 0 to 6:\n");
      printf("\n");
      fprintf_Model_Parameters(stdout, P);
      exit(0);
    }
}
