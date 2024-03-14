#include "./Include/MODEL.h"

void Label_to_Model_Parameters__LATEX(int j, char * Label, Parameter_Table *P)
{
  char * p;
  Label[0] = '\0';

  switch(j)
    {
    
    case  0:  p = strcat(Label , "Diffusion Jumping Rate");    /*  0 */
      break; 
    case  1:  p = strcat(Label , "No of INDIVIDUALS");         /*  1 */
      break;
    case  2:  p = strcat(Label , "No of CELLS");               /*  2 */
      break;
    case  3:  p = strcat(Label , "No of CELSS (X Dimension)"); /*  3 */
      break;
    case  4:  p = strcat(Label , "No of CELLS (Y Dimension)"); /*  4 */
      break;
    case  5:  p = strcat(Label , "No of RESOURCES"); /*  4 */
      break;
#ifdef DIFFUSION_HII_nD
    case  6:  p = strcat(Label, " Nu = 1/Tau\t One over the handling time (1)");  
      break;
#else
    case  6:  p = strcat(Label, "External Immigration Rate (0)");  
      break;
#endif
    case  7:  p = strcat(Label, "Decaying Rate (0)");
      break;
#ifdef DIFFUSION_HII_nD
    case  8:  p = strcat(Label, " Nu = 1/Tau\t One over the handling time (2)");  
      break;
#else
    case  8:  p = strcat(Label, "External Immigration Rate (1)");  
      break;
#endif
    case  9:  p = strcat(Label, "Decaying Rate (1)");
      break;
#ifdef DIFFUSION_AZTECA_4D
    case 10:  p = strcat(Label, "Nest Carrying Capacity (for Workers)");   /* Patch Carrying Capacity */
      break;  
#else
    case 10:  p = strcat(Label, "Patch Carrying Capacity");   /* Patch Carrying Capacity */
      break;
#endif 

    case 11: p=strcat(Label,  "Resource Local Reproduction Rate");        /* -H5 */
      break;

    case 12: p=strcat(Label,  "Consumer External Immigration Rate (0)");     /* -H5 */
      break;
    case 13: p=strcat(Label,  "Consumer Death Rate (0)");      /* -H6 */
      break;
#ifdef DIFFUSION_AZTECA_4D
    case 14: p=strcat(Label,  "Max No of Nests per Patch");     /* -H7 */
      break;
#elif DIFFUSION_AZTECA_4D_1
    case 14: p=strcat(Label,  "Jumpling/Diffusion rate of queens (1)");     /* -HuQ*/
      break;
#else
    case 14: p=strcat(Label,  "Consumer External Immigration Rate (1)");     /* -H7 */
      break;
#endif      
    case 15: p=strcat(Label,  "Consumer Death Rate (1)");      /* -H8 */ 
      break; 
      
    case 16: p=strcat(Label,  "Comsumer Attack Rate (0)");      /* -H9 */
      break;
#ifdef DIFFUSION_AZTECA_4D
    case 17: p=strcat(Label,  "Nu, Larva (Flies) Develpment Rate");  /* -H10 */
      break;
#else
    case 17: p=strcat(Label,  "Nu = 1/Tau (One over the handling time)");  /* -H10 */
      break;
#endif

    case 18: p=strcat(Label,  "Triplet formation rate (0)");        /* -H11 */
      break;
    case 19: p=strcat(Label,  "Triplet desintegration (0)");        /* -H12 */
      break;

    case 20: p=strcat(Label,  "Consumer Movement Rate (0)");        /* -H13 */
      break;

    case 21:  p = strcat(Label, "Number of Energy Levels ");   
      break;
    case 22:  p = strcat(Label, "Per Capita Fecundity: No of Offspring ");   
      break;
    case 23:  p = strcat(Label, "Energy Level at Maturity ");   
      break; 
    case 24:  p = strcat(Label, "Consummer Reproduction Rate");   
      break;
    case 25:  p = strcat(Label, "2*$k_E$ is the resourse value in energy units ");   
      break;
    case 26:  p = strcat(Label, "Energy loss maintance rate");   
      break;
    case 27:  p = strcat(Label, "Cooperation probability 1st position in the triplet");   
      break; 
    case 28:  p = strcat(Label, "Cooperation probability 2on position in the triplet");
      break;  
    case 29:  p = strcat(Label, "Establishment Rate");
      break;  
      
    default:
      printf(".... INVALID PARAMETER KEY [key = %d]\n", j);
      printf(".... The permited correspondences are:\n");
      printf("\n");
      fprintf_Model_Parameters(stdout, P);
      exit(0);
    }
}

void AssignLabel_to_Model_Parameters__LATEX(int j, char * Label, Parameter_Table * P)
{
  int i, k;
  
  char * p;
  Label[0] = '\0';
  
  if ( j < MODEL_PARAMETERS_MAXIMUM ) { 
    Label_to_Model_Parameters__LATEX(j, Label, P);
  }
  else{
    
    i = (j - MODEL_PARAMETERS_MAXIMUM)%No_of_TDC_FUNC_AUX_PARAM_MAX;
    k = (j - MODEL_PARAMETERS_MAXIMUM)/No_of_TDC_FUNC_AUX_PARAM_MAX; 
    assert( i < No_of_TDC_FUNC_AUX_PARAM_MAX && k < MODEL_PARAMETERS_MAXIMUM);
    Label_to_Model_Parameters__LATEX(k, Label, P);
    char * I_P = (char *)calloc(10, sizeof(char)); 
    sprintf(I_P, "%d", i); 
    p = strcat(Label, "(t)\\d");
    p = strcat(Label, I_P);
    p = strcat(Label, "-th Auxiliar Parameter");
    free(I_P); 
  }
  
}

void Label_to_Model_Parameters__LATEX__SYMBOL(int j, char * Label, Parameter_Table *P)
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
#ifdef DIFFUSION_HII_nD
    case  6:  p=strcat(Label,"$\\nu^{(C)}_1$");
      break;
#else
    case  6:  p = strcat(Label, "$\\lambda^{(R)}_0$");  
      break;
#endif
    case  7:  p = strcat(Label, "$\\delta^{(R)}_0$");
      break;
#ifdef DIFFUSION_HII_nD
    case  8:  p=strcat(Label,"$\\nu^{(C)}_2$");
    break;
#else
    case  8:  p = strcat(Label, "$\\lambda^{(R)}_1$");  
      break;
#endif
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
#ifdef DIFFUSION_AZTECA_4D_1
    case 14: p=strcat(Label,"$\\mu^{(Q)}$");    /* -H7 */
      break;
#else
    case 14: p=strcat(Label,"$\\lambda^{(C)}_1$");    /* -H7 */
      break;
#endif
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

    case 21:  p = strcat(Label, "$N_E$");   
      break;
    case 22:  p = strcat(Label, "$f$");   
      break;
    case 23:  p = strcat(Label, "$i_0$");   
      break; 
    case 24:  p = strcat(Label, "$\\beta_C$");   
      break;
    case 25:  p = strcat(Label, "$k_E$");   
      break;
    case 26:  p = strcat(Label, "$\\theta_C$");   
      break;
    case 27:  p = strcat(Label, "$p_1$");   
      break; 
    case 28:  p = strcat(Label, "$p_2$");   
      break;
    case 29:  p = strcat(Label, "$\\eta_R$");   
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
  int i, k; 
  char * p;
  
  Label[0] = '\0';
  if ( j < MODEL_PARAMETERS_MAXIMUM ) { 
    Label_to_Model_Parameters__LATEX__SYMBOL(j, Label, P);
  }
  else{
    
    i = (j - MODEL_PARAMETERS_MAXIMUM)%No_of_TDC_FUNC_AUX_PARAM_MAX;
    k = (j - MODEL_PARAMETERS_MAXIMUM)/No_of_TDC_FUNC_AUX_PARAM_MAX; 
    assert( i < No_of_TDC_FUNC_AUX_PARAM_MAX && k < MODEL_PARAMETERS_MAXIMUM);
    Label_to_Model_Parameters__LATEX__SYMBOL(k, Label, P);
    char * I_P = (char *)calloc(10, sizeof(char)); 
    sprintf(I_P, "%d", i); 
    p = strcat(Label, "$(t)[");
    p = strcat(Label, I_P);
    p = strcat(Label, "]$");
    free(I_P); 
  }
}
