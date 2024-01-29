#include "./Include/MODEL.h"

void Label_to_Model_Parameters(int j, char * Label, Parameter_Table *P)
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

#ifdef DIFFUSION_AZTECA_4D                                            /* -HK  */
    case 10:  p = strcat(Label, "Nest Carrying Capacity (Workers)");  /* Working Carrying Carrying Capacity per Nest */
#elif DIFFUSION_AZTECA_4D_0
    case 10:  p = strcat(Label, "Max No of Potential Nesting Trees");  /* Potential Number of Nesting Trees (per local Cell) */
#else
    case 10:  p = strcat(Label, "Patch Carrying Capacity (Resources)");   /* Patch Carrying Capacity */
      break;
#endif

    case 11: p=strcat(Label,  "Resource Local Reproduction Rate");        /* -H4 */
      break;

    case 12: p=strcat(Label,  "Consumer External Immigration Rate (0)");  /* -H5 */
      break;
    case 13: p=strcat(Label,  "Consumer Death Rate (0)");                 /* -H6 */
      break;

#ifdef DIFFUSION_AZTECA_4D
    case 14: p=strcat(Label,  "MAX No of NESTS per LOCAL PATCH");         /* -H7 */
      break;
#else      
    case 14: p=strcat(Label,  "Consumer External Immigration Rate (1)");  /* -H7 */
      break;
#endif
    case 15: p=strcat(Label,  "Consumer Death Rate (1)");                 /* -H8 */ 
      break; 
      
    case 16: p=strcat(Label,  "Comsumer Attack Rate (0)");                /* -H9  */
      break;

#ifdef DIFFUSION_AZTECA_4D
    case 17: p=strcat(Label,  "Larval Development Rate (Nu) of Flies");/* -H10 */
      break;
#elif DIFFUSION_AZTECA_4D_0 
    case 17:  p = strcat(Label, "Larval Development Rate (Nu) of Flies");  /* Working Carrying Carrying Capacity per Nest */
      break; 
#else 
    case 17: p=strcat(Label,  "Nu = 1/Tau\t One over the handling time (0)");/* -H10 */
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
    case 25:  p = strcat(Label, "2*k_E is the resourse value in energy units ");   
      break;
    case 26:  p = strcat(Label, "Maintainance Energy loss rate");   
      break;
#ifdef DIFFUSION_DRAG     
    case 27:  p = strcat(Label, "Guano Conversion Factor");   
      break;
#elif defined DIFFUSION_VGR
    case 27:  p = strcat(Label, "Guano Conversion Factor");   
      break;
#else 	      
    case 27:  p = strcat(Label, "Cooperation probability 1st position in the triplet");   
      break;
#endif      
    case 28:  p = strcat(Label, "Cooperation probability 2on position in the triplet");
      break;
    case 29:  p = strcat(Label, "Establishment Rate");
      break;
      
    default:
      printf(".... INVALID PARAMETER KEY [key = %d]\n", j);
      printf(".... The permited correspondences are (0 to 6):\n");
      printf("\n");
      fprintf_Model_Parameters(stdout, P);
      exit(0);
    }
}

void AssignLabel_to_Model_Parameters(int j, char * Label, Parameter_Table * P)
{
  int i, k;
  char * p;
  
  Label[0] = '\0';
  if ( j < MODEL_PARAMETERS_MAXIMUM ) { 
    Label_to_Model_Parameters(j, Label, P);
  }
  else{
    
    i = (j - MODEL_PARAMETERS_MAXIMUM)%No_of_TDC_FUNC_AUX_PARAM_MAX;
    k = (j - MODEL_PARAMETERS_MAXIMUM)/No_of_TDC_FUNC_AUX_PARAM_MAX; 
    assert( i < No_of_TDC_FUNC_AUX_PARAM_MAX && k < MODEL_PARAMETERS_MAXIMUM);
    Label_to_Model_Parameters(k, Label, P);
    char * I_P = (char *)calloc(10, sizeof(char)); 
    sprintf(I_P, "%d", i); 
    p = strcat(Label, " (t)[");
    p = strcat(Label, I_P);
    p = strcat(Label, "]");
    free(I_P); 
  }
}
