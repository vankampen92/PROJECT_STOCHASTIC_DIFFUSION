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
  double value;

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
    case  5: value = P->No_of_SPECIES; 
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

  return(value);
}
