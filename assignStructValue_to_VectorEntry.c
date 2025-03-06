#include <MODEL.h>

void Parameter_Table_into_Vector_Entries ( Parameter_Table * P, gsl_vector * X,
					   int * Parameter_Index, int No_of_PARAMETERS )
{
  int i;
  int key;
  double value;

  for( i=0; i<No_of_PARAMETERS; i++) {
    key = Parameter_Index[i];
    value = AssignStructValue_to_VectorEntry( key, P );
    gsl_vector_set(X, i, value);
  }
}

void Parameter_Table_into_Vector_Entries_Initial_Condition ( Parameter_Table * P, gsl_vector * X,
							     int * Index,
							     int No_of_PARAMETERS,
							     int No_of_IC)
{
  int i;
  int key;
  double value;

  for( i=0; i<No_of_IC; i++) {
    key = Index[i];
    value = Model_Variable_Initial_Condition_into_Vector_Entry_Table( key, P );
    gsl_vector_set(X, i+No_of_PARAMETERS, value);
  }
}

double AssignStructValue_to_VectorEntry(int j, Parameter_Table * P)
{
  int i, k;
  double value;

  if ( j < MODEL_PARAMETERS_MAXIMUM ) { 

    switch(j)
    {
    case  0: value = P->Mu;    
      break;
    case  1: value = P->No_of_INDIVIDUALS; 
      break;
    case  2: value = P->No_of_CELLS; 
      break;
    case  3: value = P->No_of_CELLS_X; 
      break;
    case  4: value = P->No_of_CELLS_Y; 
      break;
    case  5: value = P->No_of_RESOURCES; 
      break;
    case  6: value = P->Lambda_R_0; 
      break;
    case  7: value = P->Delta_R_0; 
      break;
    case  8: value = P->Lambda_R_1;    /* Lambda_R_1 correspond to $\gamma_0$ in LaTex: the pair formation rate*/ 
      break;
    case  9: value = P->Delta_R_1; 
      break;
    case 10: value = P->K_R;           /* Resource Carrying Capacity */ 
      break;
    case 11: value = P->Beta_R;        /* -H5 */
      break;

    case 12: value = P->Lambda_C_0;     /* -H5 */
      break;
    case 13: value = P->Delta_C_0;      /* -H6 */
      break;
    case 14: value = P->Lambda_C_1;     /* -H7 */
      break;
    case 15: value = P->Delta_C_1;      /* -H8 */ 
      break; 
      
    case 16: value = P->Alpha_C_0;      /* -H9 */
      break;
    case 17: value = P->Nu_C_0;         /* -H10 */
      break;

    case 18: value = P->Chi_C_0;        /* -H11 */
      break;
    case 19: value = P->Eta_C_0;        /* -H12 */
      break;

    case 20: value = P->Mu_C;           /* -H13 */
      break;

    case 21: value = P->N_E        ; //N_E;  /* -H14 */ /* Number of Energy Levels */
    break;
    
    case 22: value = P->f          ; //f;    /* -H15 */ /* Fecundity: Number of Offspring  */
      break;
      
    case 23: value = P->i_0        ; //i_0;  /* -H16 */ /* Energy Level at Maturity  */
      break;
      
    case 24: value = P->Beta_C     ; //Beta_C;  /* -H17 */ /* Consummer Reproduction Rate */
      break;
    
    case 25: value = P->k_E        ; /* -H18 */ /* 2* k_E is the resourse value in energy units */
      break;
      
    case 26: value = P->Theta_C    ; /* -H19 */ /* Energy loss rate for maintenance */
      break;
      
    case 27: value = P->p_1        ; /* Cooperation probability 1st position in the triplet */
      break;
      
    case 28: value = P->p_2        ; /* Cooperation probability 2on position in the triplet */  
      break;

    case 29: value = P->Eta_R     ; /* Propagule Establishment Rate  */  
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
    value = P->TDC->TDC_Function_Params[k][i]; 
  }
  
  return(value);
}
