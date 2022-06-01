#include <MODEL.h>

void Vector_Entries_into_Parameter_Table ( const gsl_vector * X, Parameter_Table * P,
					   int * Parameter_Index, int No_of_PARAMETERS )
{
  int i;
  int key;
  double value;

  for( i=0; i<No_of_PARAMETERS; i++) {
    key = Parameter_Index[i];
    value = gsl_vector_get(X, i);
    AssignVectorEntry_to_Structure( P, key, value );
  }
}

void Vector_Entries_into_Parameter_Table_Initial_Condition ( const gsl_vector * X,
							     Parameter_Table * P,
							     int * Index,
							     int No_of_PARAMETERS,
							     int No_of_IC)
{
  int i;
  int key;
  double value;

  for( i=0; i<No_of_IC; i++) {
    key = Index[i];
    value = gsl_vector_get(X, i+No_of_PARAMETERS);
    Vector_Entry_into_Model_Variable_Initial_Condition_Table( value, key, P );
  }
}

void AssignVectorEntry_to_Structure(Parameter_Table * P, int j, double value)
{
  int i,k; 
 
  if ( j < MODEL_PARAMETERS_MAXIMUM ) { 
  
    switch(j)
    {
   
    case  0: P->Mu = value;     
      break;
    case  1: P->No_of_INDIVIDUALS = value;  
      break;
    case  2: P->No_of_CELLS = value;  
      break;
    case  3: P->No_of_CELLS_X = value;  
      break;
    case  4: P->No_of_CELLS_Y = value;  
      break;
    case  5: P->No_of_RESOURCES = value;  
      break;

    case  6: P->Lambda_R_0 = value; 
      break;
    case  7: P->Delta_R_0  = value; 
      break;
    case  8: P->Lambda_R_1 = value; 
      break;
    case  9: P->Delta_R_1  = value; 
      break;
    case 10: P->K_R        = value;              /* Resource Carrying Capacity */
      break;
    case 11: P->Beta_R     = value;     /* -H4 */
      break;

    case 12: P->Lambda_C_0 = value;     /* -H5 */
      break;
    case 13: P->Delta_C_0 = value;      /* -H6 */
      break;
    case 14: P->Lambda_C_1 = value;     /* -H7 */
      break;
    case 15: P->Delta_C_1 = value;      /* -H8 */ 
      break; 
      
    case 16: P->Alpha_C_0 = value;      /* -H9 */
      break;
    case 17: P->Nu_C_0 = value;         /* -H10 */
      break;

    case 18: P->Chi_C_0 = value;        /* -H11 */
      break;
    case 19: P->Eta_C_0 = value;        /* -H12 */
      break;

    case 20: P->Mu_C = value;        /* -H12 */
      break;

    case 21: P->N_E        = value; //N_E;  /* -H14 */ /* Number of Energy Levels */
    break;
    
    case 22: P->f          = value; //f;    /* -H15 */ /* Fecundity: Number of Offspring  */
      break;
      
    case 23: P->i_0        = value; //i_0;  /* -H16 */ /* Energy Level at Maturity  */
      break;
      
    case 24: P->Beta_C     = value; //Beta_C;  /* -H17 */ /* Consummer Reproduction Rate */
      break;
      
    case 25: P->k_E        = value; /* -H18 */ /* 2* k_E is the resourse value in energy units */
      break;
      
    case 26: P->Theta_C    = value; /* -H19 */ /* Energy loss rate for maintenance */
      break;
      
    case 27: P->p_1        = value; /* Cooperation probability 1st position in the triplet */
      break;
      
    case 28: P->p_2        = value; /* Cooperation probability 2on position in the triplet */  
      break;

    case 29: P->Eta_R      = value; /* Cooperation probability 2on position in the triplet */  
      break;

    default:
      printf(".... INVALID PARAMETER KEY (key = %d)\n", j);

      printf(".... The permited correspondences are:\n");
      printf("\n");
      fprintf_Model_Parameters(stdout, P);

      printf(" The maximum number of parameters is Number_PAR_MAX\n");
      printf(" The permited number of keys go from 0, to %d\n", MODEL_PARAMETERS_MAXIMUM-1);

      exit(0);
    }
  }
  else {
    assert( P->TDC->TYPE_2_PARAMETERS > 0 && P->x_Bool == 0 );
    i = (j - MODEL_PARAMETERS_MAXIMUM)%No_of_TDC_FUNC_AUX_PARAM_MAX;
    k = (j - MODEL_PARAMETERS_MAXIMUM)/No_of_TDC_FUNC_AUX_PARAM_MAX; 
    assert( i < No_of_TDC_FUNC_AUX_PARAM_MAX && k < MODEL_PARAMETERS_MAXIMUM);
    P->TDC->TDC_Function_Params[k][i] = value; 
  } 
}
