#include "./Include/MODEL.h"

void CPGPLOT_Symbol_to_Model_Parameters(int j, char * Label, Parameter_Table *P)
{
  /* Short Labels for Model Parameters */

  char * p;
  Label[0] = '\0';

  switch(j)
    {
#ifdef DIFFUSION_ECO_1B1P
    case 0: p=strcat(Label, "$\\gm\\d0\\u");        /* -Hu */
        break;
#else 
    case  0: p=strcat(Label,"\\gm");                /* -Hu */
        break;
#endif
  
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
#ifdef DIFFUSION_HII_nD
    case  6: p=strcat(Label, "\\gn\\u(C)\\d\\d1\\u");  
      break;
#else
    case  6: p=strcat(Label, "\\gl\\u(R)\\d\\d0\\u");    
      break;
#endif
    case  7: p=strcat(Label, "\\gd\\u(R)\\d\\d0\\u");    
      break; 
#ifdef DIFFUSION_HII_nD
    case  8: p=strcat(Label, "\\gn\\u(C)\\d\\d2\\u");  
      break;
#elif DIFFUSION_ECO_1B1P
    case  8:  p = strcat(Label, "\\gg\\d0\\u"); /* $\gamma_0$ */
      break;
#else
    case  8: p=strcat(Label, "\\gl\\u(R)\\d\\d1\\u");        
      break;
#endif

    case  9: p=strcat(Label, "\\gd\\u(R)\\d\\d1\\u");    
      break; 

    case 10: p=strcat(Label, "K\\dR\\u");       
      break; 
    
    case 11: p=strcat(Label,"\\gb\\dR\\u");
      break;

    case 12: p=strcat(Label, "\\gl\\u(C)\\d\\d0\\u");    
      break;
    case 13: p=strcat(Label, "\\gd\\u(C)\\d\\d0\\u");    
      break;
#ifdef DIFFUSION_AZTECA_4D_1 
    case 14: p=strcat(Label, "\\gm\\u(Q)\\d\\d1\\u");    
      break; 
#else
    case 14: p=strcat(Label, "\\gl\\u(C)\\d\\d1\\u");    
      break;
#endif
    case 15: p=strcat(Label, "\\gd\\u(C)\\d\\d1\\u");    
      break; 
      
    case 16: p=strcat(Label, "\\ga\\u(C)\\d\\d0\\u");    
      break;

#ifdef DIFFUSION_ECO_1B1P
    case 17:  p = strcat(Label, "r");  
        break;
#else
    case 17: p=strcat(Label, "\\gn\\u(C)\\d\\d0\\u");    
      break;
#endif 

#ifdef DIFFUSION_ECO_1B1P
    case  18:  p = strcat(Label, "\\x");  /* $\chi$ in LaTeX */
      break;
#else
    case 18: p=strcat(Label, "\\gx\\u(C)\\d\\d0\\u");     
      break;
#endif

    case 19: p=strcat(Label, "\\gy\\u(C)\\d\\d0\\u");     
      break;

#ifdef DIFFUSION_ECO_1B1P
    case 20: p=strcat(Label, "$\\gm\\d1\\u");        /* -H13 */
      break;
#else 
    case 20: p=strcat(Label, "\\gm\\u(C)\\d" );      /* -H13 */
      break;
#endif
  
    case 21:  p = strcat(Label, "N\\d_E\\u");   
      break;
    case 22:  p = strcat(Label, "f");   
      break;
    case 23:  p = strcat(Label, "i\\d0\\u");   
      break; 
    case 24:  p = strcat(Label, "\\gb\\d(C)\\u");   
      break;
    case 25:  p = strcat(Label, "k\\dE\\u");   
      break;
    case 26:  p = strcat(Label, "\\gz\\d(C)\\u");   
      break;

#ifdef DIFFUSION_ECOEVO_PLANTS 
    case 27:  p = strcat(Label, "p");        /* Mutation Probability */   
      break; 
#else
    case 27:  p = strcat(Label, "p\\d1\\u");   
      break;
#endif

#ifdef DIFFUSION_ECOEVO_PLANTS
    case 28:  p = strcat(Label, "r\\d0\\u");   
      break;
#else
    case 28:  p = strcat(Label, "p\\d2\\u");   
      break;
#endif    
    
    case 29:  p = strcat(Label, "\\ge\\dR\\u");   
      break;
      
      
    default:
      printf(".... INVALID PARAMETER KEY [key=%d]\n", j);
      printf(".... The permited correspondences are (0 to 12):\n");
      printf("\n");
      fprintf_Model_Parameters(stdout, P);
      exit(0);
    }
}

void AssignCPGPLOT_Symbol_to_Model_Parameters(int j, char * Label, Parameter_Table * P)
{
  int i, k;
  char * p;
  
  Label[0] = '\0';
  if ( j < MODEL_PARAMETERS_MAXIMUM ) { 
    CPGPLOT_Symbol_to_Model_Parameters(j, Label, P);
  }
  else{
    
    i = (j - MODEL_PARAMETERS_MAXIMUM)%No_of_TDC_FUNC_AUX_PARAM_MAX;
    k = (j - MODEL_PARAMETERS_MAXIMUM)/No_of_TDC_FUNC_AUX_PARAM_MAX; 
    assert( i < No_of_TDC_FUNC_AUX_PARAM_MAX && k < MODEL_PARAMETERS_MAXIMUM);
    CPGPLOT_Symbol_to_Model_Parameters(k, Label, P);
    char * I_P = (char *)calloc(10, sizeof(char)); 
    sprintf(I_P, "%d", i); 
    p = strcat(Label, "(t)\\d");
    p = strcat(Label, I_P);
    p = strcat(Label, "\\u");
    free(I_P); 
  }
}

void Symbol_to_Model_Parameters(int j, char * Label, Parameter_Table *P)
{
  /* Short Labels for Model Parameters */

  char * p;
  Label[0] = '\0';

  switch(j)
    {

  #ifdef DIFFUSION_ECO_1B1P
      case 0: p=strcat(Label, "Mu_0");     /* -Hu */
        break;
  #else 
      case 0: p=strcat(Label, "Mu" );      /* -Hu */
        break;
  #endif
  
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

#ifdef DIFFUSION_HII_nD
    case  6: p=strcat(Label, "Nu_C_1");  
      break;
#else
    case  6: p=strcat(Label,"Lambda_R_0");    /* -H0 */
      break;
#endif

    case  7: p=strcat(Label,"Delta_R_0");     /* -H1 */  
      break; 

#ifdef DIFFUSION_HII_nD
    case  8: p=strcat(Label, "Nu_C_2");  
      break;
#elif defined DIFFUSION_ECO_1B1P
    case  8:  p = strcat(Label, "Gamma_0");   
      break;
#else
    case  8: p=strcat(Label,"Lambda_R_1");    /* -H2 */
      break;
#endif
    
    case  9: p=strcat(Label,"Delta_R_1");     /* -H3 */  
      break; 
    case 10: p=strcat(Label,"K_R");           /* -HK */
      break;
    case 11: p=strcat(Label,"Beta_R");        /* -H4 */
      break;

    case 12: p=strcat(Label,"Lambda_C_0");     /* -H5 */
      break;
    case 13: p=strcat(Label,"Delta_C_0");      /* -H6 */
      break;
#ifdef DIFFUSION_AZTECA_4D_1 
    case 14: p=strcat(Label,"Mu_Q");     /* -H7 */
      break;
#else
    case 14: p=strcat(Label,"Lambda_C_1");     /* -H7 */
      break;
#endif
    case 15: p=strcat(Label,"Delta_C_1");      /* -H8 */ 
      break; 
      
    case 16: p=strcat(Label,"Alpha_C_0");      /* -H9 */
      break;
    case 17: p=strcat(Label,"Nu_C_0");         /* -H10 */
      break;

    case 18: p=strcat(Label,"Chi_C_0");        /* -H11 */
      break;
    case 19: p=strcat(Label,"Eta_C_0");        /* -H12 */
      break;

  #ifdef DIFFUSION_ECO_1B1P
      case 20: p=strcat(Label, "Mu_1");        /* -H13 */
        break;
  #else 
      case 20: p=strcat(Label, "Mu" );      /* -H13 */
        break;
  #endif
      
    case 21:  p = strcat(Label, "N_E");   
      break;
      
    case 22:  p = strcat(Label, "f");   
      break;
      
    case 23:  p = strcat(Label, "i_0");   
      break;
      
    case 24:  p = strcat(Label, "Beta_C");   
      break;
      
    case 25:  p = strcat(Label, "k_E");   
      break;
      
    case 26:  p = strcat(Label, "Theta_C");   
      break;
 
 #ifdef DIFFUSION_ECOEVO_PLANTS     
    case 27:  p = strcat(Label, "p");   /* Mutation Probability */   
      break;
 #else
    case 27:  p = strcat(Label, "p_1");   
      break;
 #endif

#ifdef DIFFUSION_ECOEVO_PLANTS
    case 28:  p = strcat(Label, "r_0"); /* Tradeoff factor     */   
      break;
#else
    case 28:  p = strcat(Label, "p_2");   
      break;
#endif

    case 29:  p = strcat(Label, "Eta_R");   
      break;
    
    default:
      printf(".... INVALID PARAMETER KEY [key=%d]\n", j);
      printf(".... The permited correspondences are (0 to 29)\n");
      printf("\n");
      fprintf_Model_Parameters(stdout, P);
      exit(0);
    }
}

void AssignSymbol_to_Model_Parameters(int j, char * Label, Parameter_Table *P)
{
  /* Short Labels for Model Parameters */
  Label[0] = '\0';
  int i, k;
  char * p;
  
  Label[0] = '\0';
  if ( j < MODEL_PARAMETERS_MAXIMUM ) { 
    Symbol_to_Model_Parameters(j, Label, P);
  }
  else {
    
    i = (j - MODEL_PARAMETERS_MAXIMUM)%No_of_TDC_FUNC_AUX_PARAM_MAX;
    k = (j - MODEL_PARAMETERS_MAXIMUM)/No_of_TDC_FUNC_AUX_PARAM_MAX; 
    assert( i < No_of_TDC_FUNC_AUX_PARAM_MAX && k < MODEL_PARAMETERS_MAXIMUM);
    Symbol_to_Model_Parameters(k, Label, P);
    char * I_P = (char *)calloc(10, sizeof(char)); 
    sprintf(I_P, "%d", i); 
    p = strcat(Label, "(t)[");
    p = strcat(Label, I_P);
    p = strcat(Label, "]");
    free(I_P); 
  }
}
